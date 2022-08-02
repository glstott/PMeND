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
patThreshold = 0.5

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

if (params.tarball != null){
    // Input fasta files for tree building process
    input_files = Channel.fromPath( "$input_dir/*.tar" )
    log.info "List of files to be used: \n$input_files\n"
    // Process GISAID tarballs
    process untar_GISAID {
        cpus threads 
        memory "2GB"
        time "15m"
        queue "batch"
        clusterOptions "--ntasks 1"
        // This process will untar GISAID download and merge files.
        input:
        file tar from input_files

        output:
        file "${tar.simpleName}.fasta" into raw_fasta
        file "${tar.simpleName}.metadata.tsv" into raw_metadata

        script:
        """
        tar --transform "s/.*\\.metadata\\.tsv/${tar.simpleName}.metadata.tsv/" --transform "s/.*\\.fasta/${tar.simpleName}.fasta/" -xvf $tar 
        """
    }
    process collect_GISAID {
        cpus threads 
        memory "2GB"
        time "15m"
        queue "batch"
        clusterOptions "--ntasks 1"
        // This process will untar GISAID download and merge files.
        input:
        file fasta from raw_fasta.collect()
        file tsv from raw_metadata.collect()

        output:
        file "combined.fasta" into merged_fasta
        file "combined.metadata.tsv" into merged_metadata

        script:
        """
        cat *.fasta > combined.fasta
        awk 'FNR==1 && NR!=1{next;}{print}' *.metadata.tsv > combined.metadata.tsv
        """
    }

    // Split data by epiweek and generate files for IQT
    process epiweek_split {
        cpus threads 
        memory "32GB"
        time "1h"
        queue "batch"
        clusterOptions "--ntasks 1"
        cache 'lenient'

        input:
        file fasta from merged_fasta
        file "${fasta.simpleName}.metadata.tsv" from merged_metadata

        output:
        file "*.fasta" into prepped_fasta
        file "*.dt.tsv" into dates_tsv
        file "*.metadata.tsv" into metadata_tsv

        script:
        """
        python $workflow.projectDir/epiweeker.py ${fasta.simpleName}.metadata.tsv $fasta 
        """

    }
} else {
    prepped_fasta = Channel.fromPath( "$input_dir/*.fasta" ).collect()
    dates_tsv = Channel.fromPath( "$input_dir/*.tsv" ).collect()
}

// Align fasta sequences to a reference strain (Original Wuhan sequence) with MAFFT
if (run_mode != "external"){
    process mafft{
        // Initialize environment in conda
        // conda "$workflow.projectDir/envs/mafft.yaml"

        // Set slurm options.
        cpus threads 
        memory mem
        time "6h"
        queue "batch"
        clusterOptions "--ntasks $threads"
        
        // Establish output directory
        publishDir = temp_out_dir
        
        input:
        file fasta from prepped_fasta.flatten()
        
        output:
        file("${fasta.simpleName}.aligned.fasta") into alignedFasta
        
        // Add new fragments to the existing alignment set by the original wuhan sequence.
        script:
        """
        mafft --6merpair --thread ${threads} --addfragments ${fasta} $input_dir/../EPI_ISL_402124.fasta > ${fasta.simpleName}.aligned.fasta
        """

    }
}


// Check if the user is looking for fast ML trees or IQ Tree and build accordingly.
if (run_mode == 'fast'){
    process fasttree {
        // establish output directory and location for FastTree environment file
        publishDir = output_dir
        // conda "$workflow.projectDir/envs/fasttree.yaml"

        // Set slurm options. Eventually, I'll parameterize. Hard-coding for now.
        memory mem
        time "2h"
        queue "batch"
        clusterOptions "--ntasks $threads"

        input:
        file fasta from alignedFasta

        output:
        file("${fasta.name}.nwk") into phylogeny_ch
        file("${fasta.name}.nwk") into phylogeny_ch2
        
        script:
        """
        FastTree -gtr -nt $fasta > ${fasta.name}.nwk
        """
    }
} else if (run_mode == 'iqtree') {

    process iqtree {
        publishDir = output_dir
        conda "-c bioconda iqtree"
        stageInMode 'copy'

        memory mem
        time "96h" 
        queue "batch"
        clusterOptions "--ntasks $threads"

        input:
        file fasta from alignedFasta
        file "${fasta.simpleName}.dt.tsv" from dates_tsv.flatten()

        output:
        file("${fasta.simpleName}.aligned.fasta.iqt.timetree.nex") into phylogeny_ch
        file("${fasta.simpleName}.aligned.fasta.iqt.timetree.nex") into phylogeny_ch2
        
        script:
        """
        echo "${fasta.simpleName}.dt.tsv"
        iqtree2 -s $fasta -m GTR -pre ${fasta.name}.iqt -T $threads --date ${fasta.simpleName}.dt.tsv -o "hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30" -B 1000 
        """
    }

} else if (run_mode == 'external') {
    phylogeny_ch = Channel.fromPath( "$input_dir/*" )
    phylogeny_ch2 = Channel.fromPath( "$input_dir/*" )
}


// Initialize sandbox for Neo4j (using a docker image)
// https://www.nextflow.io/docs/latest/docker.html ; https://neo4j.com/developer/docker-run-neo4j/ ; https://neo4j.com/docs/operations-manual/current/docker/introduction/
// Still looking into things here. I want to ensure that the user isn't forced to start from scratch each time they begin the 
//   pipeline. It would be nice to have a NF pipeline that builds the Neo4j server for the user only if it's the first time
//   or if they don't have their own DB connection info to pass along. I'll start with the assumption of a new user each time.
// dbDir='/home/gs69042/GISAID_db'
// singularity instance start --bind /home/gs69042/conf_neo:/var/lib/neo4j/conf --bind $dbDir/data:/data --bind $dbDir/logs:/logs --bind /home/gs69042/plugins:/plugins --bind $dbDir/run:/var/lib/neo4j/run --bind $dbDir/import:/import    --env NEO4J_ACCEPT_LICENSE_AGREEMENT=yes --env NEO4J_AUTH=neo4j/test docker://neo4j:latest test_neo
// singularity shell instance://test_neo
// neo4j start
// cypher-shell -u neo4j -p test "MATCH (s:sample)-[r:calculated_distance]->(t) RETURN s,t,r LIMIT 24;"

// cypher-shell -u neo4j -p test "CREATE INDEX source_relationship3 FOR ()-[r:calculated_distance]->() ON (r.source_id);" 
// cypher-shell -u neo4j -p test "CREATE INDEX source_relationship FOR ()-[r:child_of]->() ON (r.source_id);"
requireStart = false
dbDir = '/home/gs69042/GISAID_db'
dbPath = file(dbDir)
containerSettings = "--bind /home/gs69042/conf_neo:/var/lib/neo4j/conf --bind $dbDir/data:/data --bind $dbDir/logs:/logs --bind /home/gs69042/plugins:/plugins --bind $dbDir/run:/var/lib/neo4j/run --bind $dbDir/import:/import    --env NEO4J_ACCEPT_LICENSE_AGREEMENT=yes --env NEO4J_AUTH=neo4j/test"

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
        time "30m"
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
        
        // cypher-shell -u neo4j -p test "CREATE INDEX source_relationship FOR ()-[r:child_of]->() ON (r.source_id);"
        // cypher-shell -u neo4j -p test "CREATE INDEX source_relationship3 FOR ()-[r:calculated_distance]->() ON (r.source_id);"
        script:
        """
        neo4j-admin set-initial-password test
        neo4j start
        
        until cypher-shell -u neo4j -p test "CREATE INDEX sample_name FOR (n:sample) ON (n.name);" 
        do
            echo "Create index failed, sleeping until DB running"
            sleep 6
        done
        
        cypher-shell -u neo4j -p test "CREATE INDEX phylogeny_id FOR (n:sample) ON (n.id);" 
        cypher-shell -u neo4j -p test "CREATE INDEX phylogeny_source FOR (n:LICA) ON (n.id);" 
        cypher-shell -u neo4j -p test "CREATE INDEX source_relationship FOR ()-[r:child_of]->() ON (r.source_id);"
        cypher-shell -u neo4j -p test "CREATE INDEX source_relationship3 FOR ()-[r:calculated_distance]->() ON (r.source_id);"
 
        
        echo "Mission Success"
        neo4j stop
        """
    }
} else {
    db_ch = "Existing DB"
}

// Check for user options. Patristic Only will use R. 
process tree_to_csv {
    //conda "$workflow.projectDir/envs/python.yaml"
    // change publish directory to the import folder of our database instance.
    publishDir(path: dbDir + '/import', 
        mode: 'copy', 
        overwrite: true)
    stageInMode 'copy'
    cpus threads 
    memory mem
    time "1h"
    queue "batch"
    clusterOptions "--ntasks $threads"
    // maxForks 1


    // read in the newick files and output a csv.
    input:
    file tree from phylogeny_ch

    output:
    file("${tree.simpleName}.csv") into phylo_csv

    script:
    """
    echo "Beginning python workflow step"
    python $workflow.projectDir/../data_load/treegen.py $tree ${tree.simpleName}.csv
    """
}




if (params.patristic_mode == 'R'){
    process patristic_R {
        //conda "$workflow.projectDir/envs/python.yaml"
        // change publish directory to the import folder of our database instance.
        publishDir(path: dbDir + '/import', 
            mode: 'copy', 
            overwrite: true)
        stageInMode 'copy'
        cpus threads 
        memory mem
        time "1h"
        queue "batch"
        clusterOptions "--ntasks $threads"
        // maxForks 1


        // read in the newick files and output a csv.
        input:
        file tree from phylogeny_ch2

        output:
        file("${tree.simpleName}.patristic.csv") into patristic_csv

        script:
        """
        #!/usr/bin/env Rscript
        library(ape)
        library(reshape2)
        options(digits=10)

        if (endsWith("$tree", '.nex')) tree<-read.nexus("$tree") else tree<- read.tree("$tree")
        pdm<-cophenetic(tree)
        # pp<- prop.part(tree, check.labels = TRUE)
        pdm_df <- melt(pdm,value.name='distance')
        pdm_df\$source <- "$tree"
        pdm_df <- pdm_df[pdm_df['distance'] <= as.numeric($patThreshold),]
        pdm_df <- pdm_df[pdm_df['Var1'] != pdm_df['Var2'],]
        str_cln <- function (x) {
        strsplit(strsplit(x, "\\\\|")[[1]][1], "/")[[1]][3]
        }
        pdm_df['target']<- sapply(as.character(pdm_df\$Var2), str_cln)
        pdm_df['source']<- sapply(as.character(pdm_df\$Var1), str_cln)

        write.csv(pdm_df, "${tree.simpleName}.patristic.csv", row.names=FALSE)
        """
    }

    process load_data_withPat {
        // Open up the database with docker image
        container "docker://neo4j:latest"
        containerOptions containerSettings
        
        // HPC environment settings. 
        cpus threads 
        memory mem
        time "6h"
        queue "batch"
        clusterOptions "--ntasks $threads"
        maxForks 1
        
        // Input files from phylogeny CSVs.
        input: 
        file csv from phylo_csv
        val x from db_ch
        file("${csv.simpleName}.patristic.csv") from patristic_csv

        output:
        env source_id into patristic_ch
        val "trees loaded" into tree_complete_ch

        script:
        """
        neo4j start
        
        until cypher-shell -u neo4j -p test "MATCH (n)-[r:child_of]->(m) RETURN MAX(r.source_id);" 
        do
            echo "create node failed, sleeping"
            sleep 6
        done
        
        last=\$(cypher-shell -u neo4j -p test "MATCH (n)-[r:calculated_distance]->(m) RETURN CASE WHEN MAX(r.source_id) IS NULL THEN 0 ELSE MAX(r.source_id) END;" | tail -n 1)

        if cypher-shell -u neo4j -p test "MATCH (n)-[r:calculated_distance]->(m) RETURN DISTINCT r.source" | grep -q ${csv.simpleName}; 
        then
            echo "${csv.simpleName}" "already loaded."; 
        else 
            echo "Loading new file" "${csv.simpleName}"; 

            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate( 'LOAD CSV WITH HEADERS FROM \\'file:///${csv.name}\\' AS line WITH line WHERE line.type <> \\'root\\' AND line.type <> \\'node\\' RETURN DISTINCT line.id AS id, line.type AS type', 'MERGE (child:sample:phylogeny {id: id, type: type})', {batchSize: 1000} );"
            
            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate(' LOAD CSV WITH HEADERS FROM \\'file:///${csv.simpleName}.patristic.csv\\' AS line  MATCH (child:sample:phylogeny {id: line.source})  MATCH (parent:sample {id: line.target}) WHERE id(child) < id(parent)  RETURN child, parent, line', 'CREATE (child)-[:calculated_distance {source: \\'${csv.simpleName}.patristic.csv\\', distance: toFloat(line.distance),  source_id: toInteger(\$((last+1)))}]->(parent)', {batchSize: 1000});"
            
        fi
        source_id=\$((last+1))
        echo "Mission Success"
        neo4j stop
        """
    }

    
} else {
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
        maxForks 1
        
        // Input files from phylogeny CSVs.
        input: 
        file csv from phylo_csv
        val x from db_ch

        output:
        env source_id into loaded_tree_ch
        val "trees loaded" into tree_complete_ch

        script:
        """
        neo4j start
        
        until cypher-shell -u neo4j -p test "MATCH (n)-[r:child_of]->(m) RETURN MAX(r.source_id);" 
        do
            echo "create node failed, sleeping"
            sleep 6
        done
        
        last=\$(cypher-shell -u neo4j -p test "MATCH (n)-[r:child_of]->(m) RETURN CASE WHEN MAX(r.source_id) IS NULL THEN 0 ELSE MAX(r.source_id) END;" | tail -n 1)

        if cypher-shell -u neo4j -p test "MATCH (n)-[r:child_of]->(m) RETURN DISTINCT r.source" | grep -q ${csv.simpleName}; 
        then
            echo "${csv.simpleName}" "already loaded."; 
        else 
            echo "Loading new file" "${csv.simpleName}"; 

            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate('LOAD CSV WITH HEADERS FROM \\'file:///${csv.simpleName}\\' AS line WITH line WHERE line.type = \\'root\\' OR line.type = \\'node\\' RETURN DISTINCT line.id AS id, line.type AS type, line.longid AS longid, line.source AS source', 'MERGE (child:LICA:phylogeny {id: id, long_id: longid, type: type})',     {batchSize: 1000} );"

            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate( 'LOAD CSV WITH HEADERS FROM \\'file:///${csv.simpleName}\\' AS line WITH line WHERE line.type <> \\'root\\' AND line.type <> \\'node\\' RETURN DISTINCT line.id AS id, line.type AS type', 'MERGE (child:sample:phylogeny {id: id, type: type})', {batchSize: 1000} );"

            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate(' LOAD CSV WITH HEADERS FROM \\'file:///${csv.simpleName}\\' AS line  WITH line WHERE line.type = \\'leaf\\'  MATCH (child:sample:phylogeny {id: line.id})  MATCH (parent:LICA:phylogeny {id: line.parent})  RETURN child, parent, line', 'CREATE (child)-[:child_of {source: \\'${csv.simpleName}\\', distance: toFloat(line.length),  source_id: toInteger(\$((last+1)))}]->(parent)', {batchSize: 1000});"

            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate('LOAD CSV WITH HEADERS FROM \\'file:///${csv.simpleName}\\' AS line  WITH line  WHERE line.type <> \\'leaf\\' AND line.type <> \\'root\\'  MATCH (child:LICA:phylogeny {id: line.id})  MATCH (parent:LICA:phylogeny {id: line.parent})  RETURN child, parent, line ', 'CREATE (child)-[:child_of {source: \\'${csv.simpleName}\\', distance: toFloat(line.length),  source_id:  toInteger(\$((last+1)))}]->(parent) ',  {batchSize: 1000}); "
            
        fi
        source_id=\$((last+1))
        echo "Mission Success"
        neo4j stop
        """
    }
    
    process patristic_calculation {
        // Spawn docker image, accepting license agreement and binding relevant paths. 
        container "docker://neo4j:latest"
        containerOptions containerSettings

        // Specify parameters for our HPC environment. We should move these to a conf file in the future to clean things up. 
        cpus threads 
        memory mem
        time "6h"
        queue "batch"
        clusterOptions "--ntasks $threads"
        maxForks 1

        input:
        val source_id from loaded_tree_ch.collect().flatten()

        output:
        val "$source_id" into patristic_ch

        script:
        """
        neo4j start
        
        until cypher-shell -u neo4j -p test "MATCH (n)-[r:child_of]->(m) RETURN DISTINCT r.source_id" 
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

            echo "CALL gds.graph.project('trees', 'phylogeny', {child_of: {properties: ['source_id', 'distance'], orientation: 'UNDIRECTED'}}) YIELD graphName AS gn, nodeCount AS nc, relationshipCount AS rc WITH gn, nc, rc CALL gds.beta.graph.project.subgraph('one_tree', 'trees', '*', 'r.source_id > ${source_id.toInteger()-1}.0 AND r.source_id < ${source_id.toInteger()+1}.0') yield graphName AS subgraph_name, nodeCount AS subgraph_nc, relationshipCount AS subgraph_rc RETURN gn, nc, rc, subgraph_name, subgraph_nc, subgraph_rc; "

            cypher-shell -u neo4j -p test "CALL gds.graph.project('trees', 'phylogeny', {child_of: {properties: ['source_id', 'distance'], orientation: 'UNDIRECTED'}}) YIELD graphName AS gn, nodeCount AS nc, relationshipCount AS rc WITH gn, nc, rc CALL gds.beta.graph.project.subgraph('one_tree', 'trees', '*', 'r.source_id > ${source_id.toInteger()-1}.0 AND r.source_id < ${source_id.toInteger()+1}.0') yield graphName AS subgraph_name, nodeCount AS subgraph_nc, relationshipCount AS subgraph_rc RETURN gn, nc, rc, subgraph_name, subgraph_nc, subgraph_rc; "

            echo "CALL apoc.periodic.iterate('MATCH (source:sample), (target:sample) WITH source, target WHERE id(source) < id(target) AND (source)-[:child_of {source_id: $source_id}]->() AND (target)-[:child_of {source_id: $source_id}]->() RETURN source, target', ' CALL gds.shortestPath.dijkstra.stream(\\'one_tree\\', {  relationshipWeightProperty: \\'distance\\', sourceNode: id(source), targetNode: id(target) }) YIELD sourceNode, targetNode, totalCost AS distance WHERE gds.util.isFinite(distance) = true AND distance <= $patThreshold WITH source, target, distance CREATE (source)-[:calculated_distance {distance: distance,  source_id: $source_id}]->(target)', {batchSize: 1000});"
            
            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate('MATCH (source:sample), (target:sample) WITH source, target WHERE id(source) < id(target) AND (source)-[:child_of {source_id: $source_id}]->() AND (target)-[:child_of {source_id: $source_id}]->() RETURN source, target', ' CALL gds.shortestPath.dijkstra.stream(\\'one_tree\\', {  relationshipWeightProperty: \\'distance\\', sourceNode: id(source), targetNode: id(target) }) YIELD sourceNode, targetNode, totalCost AS distance WHERE gds.util.isFinite(distance) = true AND distance <= $patThreshold WITH source, target, distance CREATE (source)-[:calculated_distance {distance: distance,  source_id: $source_id}]->(target)', {batchSize: 1000});"
            
        fi
        
        echo "Mission Success"
        neo4j stop
        """

    }
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
        time "48h"
        queue "batch"
        clusterOptions "--ntasks $threads"
        maxForks 1

        input:
        val source_id from patristic_ch.collect().flatten()

        script:
        """
        neo4j start
        
        until cypher-shell -u neo4j -p test "MATCH (n)-[r]->(m) RETURN MAX(r.source_id);"
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

            cypher-shell -u neo4j -p test "MATCH (n:sample)-[:calculated_distance {source_id: $source_id}]->()  WITH DISTINCT id(n) AS iteration_number LIMIT $MST_number CALL gds.alpha.spanningTree.minimum.write('subgraph', {startNodeId: iteration_number,    relationshipWeightProperty: 'distance', writeProperty: 'MST_' + iteration_number + '_' + '$source_id', weightWriteProperty: 'distance_combined'}) YIELD computeMillis, writeMillis RETURN computeMillis AS coM, writeMillis AS wM;"
            
            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate(    'MATCH (n:sample)-[r]-(:sample) WHERE type(r) ENDS WITH \\'$source_id\\' WITH n, COUNT(DISTINCT type(r)) AS max_width  MATCH (n)-[r]->(n2:sample) WHERE type(r) ENDS WITH \\'$source_id\\'  AND id(n) < id(n2)  RETURN max_width, n, n2, COUNT(DISTINCT r) AS edge_width, AVG(r.distance_combined) AS mean_dist', 'CREATE (n)-[:joint_MST {source: \\'$source_id\\', edge_width: edge_width / toFloat(max_width), mean_dist: mean_dist}]->(n2)', {batchSize: 10000}) YIELD batches AS batches, timeTaken AS timeTaken, failedBatches AS failedBatches RETURN batches, timeTaken, failedBatches;"

            cypher-shell -u neo4j -p test "CALL apoc.periodic.iterate('MATCH (:sample)-[r]->(:sample) WHERE type(r) ENDS WITH \\'$source_id\\' RETURN r', 'DELETE r', {batchSize: 10000}) YIELD batches AS batches, timeTaken AS timeTaken, failedBatches AS failedBatches RETURN batches, timeTaken, failedBatches;"

        fi
        
        echo "Mission Success"
        neo4j stop
        """

    }
}