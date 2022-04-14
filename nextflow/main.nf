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
        cpus threads
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

        cpus threads
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