#!/usr/bin/env nextflow
// Authors: Garrick Stott (project leader), Noah Legall 
// Purpose: A simple use case of using Nextflow to populate a Tree Aligned Graph (TAG)

// Print log info to user screen.
log.info """ 
    PMeND (v0.1)    
=============================
A graph database which enables storage of phylogenies as 
Tree Aligned Graphs (TAGs) and integrates these data with sample metadata. 

Project : $workflow.projectDir
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
Manifest's pipeline version: $workflow.manifest.version
=============================
"""

// Set default parameter values
temp_out_dir = "./"
output_dir = "../out/"
mem = "64GB"
threads = 8
run_mode = "fast"
input_dir = "/scratch/gs69042/PMeND/data/example"

// Check for user inputs
if (params.input != null){
   input_dir = params.input 
} 
if (params.temp_out_dir != null){
    temp_out_dir = params.temp_out_dir
}
if (params.output_dir != null){
    output_dir = params.output_dir
}
if (params.run_mode != null){
    run_mode = params.run_mode
}
threads = 1
if (params.threads != null){
    threads = params.threads
}

// Input fasta files for tree building process
input_files = Channel.fromPath( "$input_dir/*.fasta" )
log.info "List of files to be used: \n$input_files\n"

// Align fasta sequences to a reference strain (Original Wuhan sequence) with MAFFT
process mafft{
    // Initialize environment in conda
    conda "$workflow.projectDir/envs/mafft.yaml"

    // Set slurm options.
    cpus threads 
    memory mem
    time "6h"
    queue "batch"
    clusterOptions "--ntasks $threads"
    
    // Establish output directory
    publishDir = temp_out_dir
    
    input:
    file fasta from input_files
    
    output:
    file("${fasta.simpleName}.aligned.fasta") into alignedFasta
    
    // Add new fragments to the existing alignment set by the original wuhan sequence.
    script:
    """
    mafft --6merpair --thread ${threads} --addfragments ${fasta} $input_dir/../EPI_ISL_402124.fasta > ${fasta.simpleName}.aligned.fasta
    """

}


// Check if the user is looking for fast ML trees or IQ Tree and build accordingly.
if (run_mode == 'fast'){
    process fasttree {
        // establish output directory and location for FastTree environment file
        publishDir = output_dir
        conda "$workflow.projectDir/envs/fasttree.yaml"

        // Set slurm options. Eventually, I'll parameterize. Hard-coding for now.
        memory mem
        time "6h"
        queue "batch"
        clusterOptions "--ntasks $threads"

        input:
        file fasta from alignedFasta

        output:
        file("${fasta.name}.nwk") into phylogeny_ch
        
        script:
        """
        FastTree -gtr -nt $fasta > ${fasta.name}.nwk
        """
    }
} else {

    process iqtree {
        publishDir = output_dir
        conda "$workflow.projectDir/envs/iqtree.yaml"

        memory mem
        time "72h"
        queue "batch"
        clusterOptions "--ntasks $threads"

        input:
        file fasta from alignedFasta

        output:
        file("*") into phylogeny_ch
        
        script:
        """
        iqtree -s $fasta -m GTR -pre ${fasta.name}.iqt -nt $threads
        """
    }

}


// Initialize sandbox for Neo4j (using a docker image)
// https://www.nextflow.io/docs/latest/docker.html ; https://neo4j.com/developer/docker-run-neo4j/ ; https://neo4j.com/docs/operations-manual/current/docker/introduction/
// Still looking into things here. I want to ensure that the user isn't forced to start from scratch each time they begin the 
//   pipeline. It would be nice to have a NF pipeline that builds the Neo4j server for the user only if it's the first time
//   or if they don't have their own DB connection info to pass along. I'll start with the assumption of a new user each time.
// singularity instance start     --bind $HOME/neo4j_empty/data:/data     --bind $HOME/neo4j_empty/logs:/logs     --bind $HOME/neo4j_empty/import:/import     --bind $HOME/plugins:/plugins --bind /home/gs69042/conf_neo:/var/lib/neo4j/conf --bind /home/gs69042/neo4j_empty/run:/var/lib/neo4j/run   --env NEO4J_AUTH=neo4j/test --env NEO4J_ACCEPT_LICENSE_AGREEMENT=yes docker://neo4j:latest test_neo
// singularity shell instance://test_neo
// neo4j start
// cypher-shell -u neo4j -p test "MATCH (s:sample)-[r:calculated_distance]->(t) RETURN s,t,r LIMIT 24;"

// cypher-shell -u neo4j -p test "CREATE INDEX source_relationship3 FOR ()-[r:calculated_distance]->() ON (r.source_id);" 
// cypher-shell -u neo4j -p test "CREATE INDEX source_relationship FOR ()-[r:child_of]->() ON (r.source_id);"
requireStart = false
dbDir = '/home/gs69042/neo4j_empty'
dbPath = file(dbDir)
containerSettings = "--bind /home/gs69042/conf_neo:/var/lib/neo4j/conf --bind /home/gs69042/neo4j_empty/data:/data --bind /home/gs69042/neo4j_empty/logs:/logs --bind /home/gs69042/plugins:/plugins --bind /home/gs69042/neo4j_empty/run:/var/lib/neo4j/run --bind /home/gs69042/neo4j_empty/import:/import    --env NEO4J_ACCEPT_LICENSE_AGREEMENT=yes --env NEO4J_AUTH=neo4j/test"

// If we don't see folders in the target directory, go ahead and make them. 
//   Current state is naive, checking for quantity rather than specific paths. I will come back and correct @ a later date.
if (dbPath.list().size() == 0) {
    // Create subdirectories needed for Neo4j Docker Image. Singularity doesn't allow for auto-create, and pre-scripts don't allow for containerless script generation.  
    requireStart = true 
    process create_db_env {
        cpus threads 
        memory "1GB"
        time "15m"
        queue "batch"
        clusterOptions "--ntasks 1"
        // we mustn't cache this one. Otherwise if something happens to the DB offline, it won't be correctly reinstantiated. This runs parallel with other processes so it shouldn't be a problem.
        cache false

        output:
        val "Environment Ready" into sequential_ch

        shell:
        """
        mkdir $dbDir/logs
        mkdir $dbDir/run
        mkdir $dbDir/import
        mkdir $dbDir/data
        """
    }

    process new_db {
        // Spawn docker image, accepting license agreement and binding relevant paths. 
        container "docker://neo4j:latest"
        containerOptions containerSettings

        // Specify parameters for our HPC environment. We should move these to a conf file in the future to clean things up. 
        cpus threads 
        memory mem
        time "15m"
        queue "batch"
        clusterOptions "--ntasks $threads"
        // we mustn't cache this one. Otherwise if something happens to the DB offline, it won't be correctly reinstantiated. This runs parallel with other processes so it shouldn't be a problem.
        cache false

        input:
        val x from sequential_ch

        output:
        val "Step complete" into db_ch
        
        // To spin up the database for the first time we:
        //   1. set initial password
        //   2. Spin up database and generate indices where needed downstream
        script:
        """
        neo4j-admin set-initial-password test
        neo4j start
        
        until cypher-shell -u neo4j -p test "CREATE INDEX sample_name FOR (n:sample) ON (n.name);" 
        do
            echo "Create index failed, sleeping until DB running"
            sleep 6
        done
        
        cypher-shell -u neo4j -p test "CREATE INDEX phylogeny_id FOR (n:phylogeny) ON (n.id);" 
 
        
        echo "Mission Success"
        neo4j stop
        """
    }
}

process tree_to_csv {
    //conda "$workflow.projectDir/envs/python.yaml"
    // change publish directory to the import folder of our database instance.
    publishDir(path: dbDir + '/import', 
        mode: 'copy', 
        overwrite: true)
    memory mem
    time "1h"
    queue "batch"
    clusterOptions "--ntasks $threads"

    // read in the newick files and output a csv.
    input:
    file phylogeny from phylogeny_ch

    output:
    file("${phylogeny.name}.csv") into phylo_csv

    script:
    """
    python $workflow.projectDir/../data_load/treegen.py $phylogeny ${phylogeny.name}.csv
    """
}

process load_data {
    // Open up the database with docker image
    container "docker://neo4j:latest"
    containerOptions containerSettings
    
    // HPC environment settings. 
    cpus threads 
    memory mem
    time "6h"
    queue "batch"
    clusterOptions "--ntasks $threads"
    
    // Input files from phylogeny CSVs.
    input: 
    file csv from phylo_csv
    if (requireStart) {
        val entry from db_ch
    }

    output:
    env source_id into loaded_tree_ch

    script:
    """
    neo4j start
    
    until cypher-shell -u neo4j -p test "MATCH (n)-[r:child_of]->(m) RETURN MAX(r.source_id);" 
    do
        echo "create node failed, sleeping"
        sleep 6
    done
    
    last=\$(cypher-shell -u neo4j -p test "MATCH (n)-[r:child_of]->(m) RETURN MAX(r.source_id);" | tail -n 1)

    if cypher-shell -u neo4j -p test "MATCH (n)-[r:child_of]->(m) RETURN DISTINCT r.source" | grep -q ${csv.name}; 
    then
        echo "${csv.name}" "already loaded."; 
    else 
        echo "Loading new file" "${csv.name}"; 

        cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate('LOAD CSV WITH HEADERS FROM \\'file:///${csv.name}\\' AS line WITH line WHERE line.type = \\'root\\' OR line.type = \\'node\\' RETURN DISTINCT line.id AS id, line.type AS type', 'MERGE (child:LICA:phylogeny {id: id, type: type})',     {batchSize: 1000} );"

        cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate( 'LOAD CSV WITH HEADERS FROM \\'file:///${csv.name}\\' AS line WITH line WHERE line.type <> \\'root\\' AND line.type <> \\'node\\' RETURN DISTINCT line.id AS id, line.type AS type', 'MERGE (child:sample:phylogeny {id: id, type: type})', {batchSize: 1000} );"

        cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate(' LOAD CSV WITH HEADERS FROM \\'file:///${csv.name}\\' AS line  WITH line WHERE line.type = \\'leaf\\'  MATCH (child:sample:phylogeny {id: line.id})  MATCH (parent:LICA:phylogeny {id: line.parent})  RETURN child, parent, line', 'CREATE (child)-[:child_of {source: \\'${csv.name}\\', distance: toFloat(line.length),  source_id: toInteger(\${last+1})}]->(parent)', {batchSize: 1000});"

        cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate('LOAD CSV WITH HEADERS FROM \\'file:///${csv.name}\\' AS line  WITH line  WHERE line.type <> \\'leaf\\' AND line.type <> \\'root\\'  MATCH (child:LICA:phylogeny {id: line.id})  MATCH (parent:LICA:phylogeny {id: line.parent})  RETURN child, parent, line ', 'CREATE (child)-[:child_of {source: \\'${csv.name}\\', distance: toFloat(line.length),  source_id:  toInteger(\${last+1})}]->(parent) ',  {batchSize: 1000}); "
        
    fi
    source_id=\${last+1}
    echo "Mission Success"
    neo4j stop
    """
}

process patristicCalculation {
    // Spawn docker image, accepting license agreement and binding relevant paths. 
    container "docker://neo4j:latest"
    containerOptions containerSettings

    // Specify parameters for our HPC environment. We should move these to a conf file in the future to clean things up. 
    cpus threads 
    memory mem
    time "15m"
    queue "batch"
    clusterOptions "--ntasks $threads"

    input:
    val source_id from loaded_tree_ch

    output:
    val "$source_id" into patristic_ch

    script:
    """
    neo4j start
    
    until cypher-shell -u neo4j -p test "MATCH (n)-[r:child_of]->(m) RETURN MAX(r.source_id);"
    do
        echo "create node failed, sleeping"
        sleep 6
    done
    
    if cypher-shell -u neo4j -p test "MATCH (n)-[r:calculated_distance]->(m) RETURN DISTINCT r.source_id" | grep -q ${source_id}; 
    then
        echo "${source_id}" "already loaded."; 
    else 
        echo "Loading new patristic distances" "${source_id}"; 

        cypher-shell -u neo4j -p test "CALL gds.graph.drop('trees', false);"
        cypher-shell -u neo4j -p test "CALL gds.graph.drop('one_tree', false);"

        cypher-shell -u neo4j -p test "CALL gds.graph.project('trees', 'phylogeny', {child_of: {properties: ['source_id', 'distance'], orientation: 'UNDIRECTED'}}) YIELD graphName AS gn, nodeCount AS nc, relationshipCount AS rc WITH gn, nc, rc CALL gds.beta.graph.project.subgraph('one_tree', 'trees', '*', 'r.source_id > ${source_id.toInteger()-1}.0 AND r.source_id < ${source_id.toInteger()+1}.0') yield graphName AS subgraph_name, nodeCount AS subgraph_nc, relationshipCount AS subgraph_rc RETURN gn, nc, rc, subgraph_name, subgraph_nc, subgraph_rc; "

        cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate('MATCH (source:sample)-[r:child_of {source_id: $source_id}]->(), (target:sample)-[r2:child_of {source_id: $source_id}]->() WITH source, target, r.source AS src WHERE id(source) < id(target) AND r.source=r2.source CALL gds.shortestPath.dijkstra.stream(\\'one_tree\\', {  relationshipWeightProperty: \\'distance\\', sourceNode: id(source), targetNode: id(target) }) YIELD sourceNode, targetNode, totalCost AS distance WHERE gds.util.isFinite(distance) = true AND distance <= 0.5 RETURN source, target, distance, src', 'CREATE (source)-[:calculated_distance {source: src, distance: distance,  source_id: $source_id}]->(target)', {batchSize: 1000});"
        
    fi
    
    echo "Mission Success"
    neo4j stop
    """

}

MST_number=2000
if ((params.mst != null) && (params.mst)) {
    process generateMST {
        // Spawn docker image, accepting license agreement and binding relevant paths. 
        container "docker://neo4j:latest"
        containerOptions containerSettings

        // Specify parameters for our HPC environment. We should move these to a conf file in the future to clean things up. 
        cpus threads 
        memory mem
        time "15m"
        queue "batch"
        clusterOptions "--ntasks $threads"

        input:
        val source_id from patristic_ch

        script:
        """
        neo4j start
        
        until cypher-shell -u neo4j -p test "MATCH (n)-[r:child_of]->(m) RETURN MAX(r.source_id);"
        do
            echo "create node failed, sleeping"
            sleep 6
        done
        
        if cypher-shell -u neo4j -p test "MATCH (n)-[r:joint_MST]->(m) RETURN DISTINCT r.source_id" | grep -q ${source_id}; 
        then
            echo "${source_id}" "already loaded."; 
        else 
            echo "Loading new patristic distances" "${source_id}"; 

            cypher-shell -u neo4j -p test "CALL gds.graph.drop('distance_graph', false);"
            cypher-shell -u neo4j -p test "CALL gds.graph.drop('subgraph', false);"

            cypher-shell -u neo4j -p test "CALL gds.graph.project('distance_graph', 'sample', {calculated_distance: {orientation: 'UNDIRECTED', properties: ['source_id', 'distance']}} ) YIELD graphName AS gn, nodeCount AS nc, relationshipCount AS rc RETURN gn, nc, rc;"

            cypher-shell -u neo4j -p test "CALL gds.beta.graph.project.subgraph('subgraph', 'distance_graph', '*', 'r.source_id > ${source_id.toInteger()-1}.0 AND r.source_id < ${source_id.toInteger()+1}.0'  ) YIELD graphName AS gn, nodeCount AS nc, relationshipCount AS rc RETURN gn, nc, rc;" 

            cypher-shell -u neo4j -p test "MATCH (n:sample)-[:child_of {source_id: $source_id}]->()  WITH DISTINCT id(n) AS iteration_number LIMIT $MST_number CALL gds.alpha.spanningTree.minimum.write('subgraph', {startNodeId: iteration_number,    relationshipWeightProperty: 'distance', writeProperty: 'MST_' + iteration_number + '_' + '$source_id', weightWriteProperty: 'distance_combined'}) YIELD computeMillis, writeMillis RETURN computeMillis AS coM, writeMillis AS wM;"
            
            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate(    'MATCH (n:sample)-[r]-(:sample) WHERE type(r) ENDS WITH \\'$source_id\\' WITH n, COUNT(DISTINCT type(r)) AS max_width  MATCH (n)-[r]->(n2:sample) WHERE type(r) ENDS WITH \\'$source_id\\'  AND id(n) < id(n2)  RETURN max_width, n, n2, COUNT(DISTINCT r) AS edge_width, AVG(r.distance_combined) AS mean_dist', 'CREATE (n)-[:joint_MST {source: \\'$source_id\\', edge_width: edge_width / toFloat(max_width), mean_dist: mean_dist}]->(n2)', {batchSize: 10000}) YIELD batches AS batches, timeTaken AS timeTaken, failedBatches AS failedBatches RETURN batches, timeTaken, failedBatches;"

            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate('MATCH (:sample)-[r]->(:sample) WHERE type(r) ENDS WITH \\'$source_id\\' RETURN r', 'DELETE r', {batchSize: 10000}) YIELD batches AS batches, timeTaken AS timeTaken, failedBatches AS failedBatches RETURN batches, timeTaken, failedBatches;"

        fi
        
        echo "Mission Success"
        neo4j stop
        """

    }
}