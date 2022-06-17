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
run_mode = "iqtree"

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
input_files = Channel.fromPath( "$input_dir*.fasta" )
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
        file("${fasta.baseName}.nwk") into phylogeny_ch
        
        script:
        """
        FastTree -gtr -nt $fasta > ${fasta.baseName}.nwk
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
        iqtree -s $fasta -p /scratch/gs69042/PMeND/phylogeny_generation/partitions.txt -m GTR -pre ${fasta.baseName}.iqt -nt $threads
        """
    }

}


// Initialize sandbox for Neo4j (using a docker image)
// https://www.nextflow.io/docs/latest/docker.html ; https://neo4j.com/developer/docker-run-neo4j/ ; https://neo4j.com/docs/operations-manual/current/docker/introduction/
// Still looking into things here. I want to ensure that the user isn't forced to start from scratch each time they begin the 
//   pipeline. It would be nice to have a NF pipeline that builds the Neo4j server for the user only if it's the first time
//   or if they don't have their own DB connection info to pass along. I'll start with the assumption of a new user each time.
// singularity instance start     --bind $HOME/neo4j_empty/data:/data     --bind $HOME/neo4j_empty/logs:/logs     --bind $HOME/neo4j_empty/import:/var/lib/neo4j/import     --bind $HOME/neo4j_empty/plugins:/plugins --bind $HOME/neo4j_empty/conf:/var/lib/neo4j/conf    --env NEO4J_AUTH=neo4j/test --env NEO4J_ACCEPT_LICENSE_AGREEMENT=yes docker://neo4j:latest test_neo
// singularity shell instance://test_neo
// neo4j start
// cypher-shell -u neo4j -p test
dbDir = '/home/gs69042/neo4j_empty'
dbPath = file(dbDir)
if (dbPath.list().size() == 0) {
    // Create subdirectories needed for Neo4j Docker Image
    new File(dbDir + '/data').mkdir()
    new File(dbDir + '/logs').mkdir()
    new File(dbDir + '/run').mkdir()
    new File(dbDir + '/plugins').mkdir()
    new File(dbDir + '/import').mkdir()
    
    process new_db {
        container "docker://neo4j:latest"
        containerOptions "--bind /home/gs69042/conf_neo:/var/lib/neo4j/conf --bind /home/gs69042/neo4j_empty/data:/data --bind /home/gs69042/neo4j_empty/logs:/logs --bind /home/gs69042/neo4j_empty/plugins:/plugins --bind /home/gs69042/neo4j_empty/run:/var/lib/neo4j/run    --env NEO4J_ACCEPT_LICENSE_AGREEMENT=yes "

        cpus threads 
        memory mem
        time "6h"
        queue "batch"
        clusterOptions "--ntasks $threads"
        
        script:
        """
        neo4j-admin set-initial-password test

        curl https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/4.3.0.6/apoc-4.3.0.6-all.jar > $dbDir/plugins/apoc-4.3.0.6-all.jar
        curl https://github.com/neo4j/graph-data-science/releases/download/2.1.2/neo4j-graph-data-science-2.1.2.jar > $dbDir/plugins/neo4j-graph-data-science-2.1.2.jar

        neo4j start
        
        until cypher-shell -u neo4j -p test "CREATE INDEX sample_name FOR (n:sample) ON (n.name);" 
        do
            echo "create node failed, sleeping"
            sleep 6
        done
        
        cypher-shell -u neo4j -p test "CREATE INDEX phylogeny_id FOR (n:phylogeny) ON (n.id);" 
        cypher-shell -u neo4j -p test "CREATE INDEX source_relationship3 FOR ()-[r:calculated_distance]->() ON (r.source_id);" 
        cypher-shell -u neo4j -p test "CREATE INDEX source_relationship FOR ()-[r:child_of]->() ON (r.source_id);" 
        
        echo "Mission Success"
        neo4j stop
        """
    }
}

process tree_to_csv {
    publishDir = dbDir + '/import'
    input:
    file phylogeny from phylogeny_ch
    output:
    file("${phylogeny.baseName}.csv") into phylo_csv

    script:
    '''
    python3 ../data_load/treegen.py $phylogeny ${phylogeny.baseName}.csv
    '''
}

process load_data {
    container "docker://neo4j:latest"
    containerOptions "--bind /home/gs69042/conf_neo:/var/lib/neo4j/conf --bind /home/gs69042/neo4j_empty/data:/data --bind /home/gs69042/neo4j_empty/logs:/logs --bind /home/gs69042/neo4j_empty/plugins:/plugins --bind /home/gs69042/neo4j_empty/run:/var/lib/neo4j/run    --env NEO4J_ACCEPT_LICENSE_AGREEMENT=yes "
    

    cpus threads 
    memory mem
    time "6h"
    queue "batch"
    clusterOptions "--ntasks $threads"
    
    input: 
    file csv from phylo_csv

    script:
    """
    neo4j start
    
    until cypher-shell -u neo4j -p test "CREATE (t:test {id: 'AbsolutelyRidiculous'});" 
    do
        echo "create node failed, sleeping"
        sleep 6
    done

    
    echo "Mission Success"
    neo4j stop
    """
}

process analyze_tree {

}