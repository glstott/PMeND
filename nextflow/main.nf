#!/usr/bin/env nextflow
// Authors: Garrick Stott (project leader), Noah Legall 
// Purpose: A simple use case of using Nextflow to populate a Tree Aligned Graph (TAG)
//          Potential to automate fully 

// We start off with a message that prints out to the screen
// Useful from a user's perspective - at least we see basic info on the pipeline
// the words in $'s are variables, and the ones referenced here are related to the nextflow script in general
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

// So we start by having to establish where the data comes from. 
// As a default location we can assume maybe it's the current directory
temp_out_dir = "./"
output_dir = "../out/"
mem = "64GB"

// but we can add flexibility by allowing the user provide a parameter
// Of course, we have to check if the parameter variable is populated
if (params.input != null){
   input_dir = params.input 
}

if (params.temp_out_dir != null){
    temp_out_dir = params.temp_out_dir
}

// As for other parameters, such as number of threads, we can make that a user
// defined parameter too. Use what you did above and define a thread parameter
threads = 1
if (params.threads != null){
    threads = params.threads
}

// Nextflow has two really important features - the 'process' and 'channels'
// Channels are First in, First out datastructures that process individual data in 
// specific ways. Here, we use a simple regex to find the particular files that will be 
// inputs for your process
input_files = Channel.fromPath( "$input_dir*.fasta" )

log.info "$input_files"

process mafft{
    
    // A cool thing about Nextflow is that we can treat each process like it's  own little environment
    // Conda is a tool that allows use to specify computing environments without leading to dependency clashes. 
    conda "$workflow.projectDir/envs/mafft.yaml"//where can we get the information on the environment from? maybe check out mbovpan for an idea

    cpus threads //perfect place to put your thread parameter
    memory mem
    time "6h"
    queue "batch"
    clusterOptions "--ntasks 16"
    
    publishDir = temp_out_dir

    // There might be other important process attributes to add, but I'll leave that to you to discover -
    // https://www.nextflow.io/docs/latest/process.html
    
    // The process needs files to work on, so we use the input directive to establish the files to work on
    // it pulls from our initial channel.
    input:
    file fasta from input_files

    // Here we need to establish an output that will be funneled into a new channel
    output:
    file("${fasta.simpleName}.aligned.fasta") into alignedFasta
    // from your command, what is output that we need to work with in future processes?
    // look at the script portion below, maybe the output is specified there.

    // Last but not least, we need to give instructions on what we want to run!
    // below you need to change a few things to have it incorporate the data the process specifies
    // for example, how do we incorporate the thread parameter? or the input files? Also where are we getting the 'reference.fasta'
    // give it a try to figure it out!
    script:
    """
    mafft --6merpair --thread ${threads} --addfragments ${fasta} /scratch/gs69042/PMeND/data/EPI_ISL_402124.fasta > ${fasta.simpleName}.aligned.fasta
    """

}

// And that's kind of the basics! from there, you can solutions by googling
// heres a more blank version of the iqtree process after aligning the sequences

process fasttree {
    publishDir = output_dir

    conda "$workflow.projectDir/envs/fasttree.yaml"

    cpus threads
    memory mem
    time "6h"
    queue "batch"
    clusterOptions "--ntasks 16"



    input:
    file fasta from alignedFasta

    output:
    file("${fasta.baseName}.nwk")
    
    script:
    """
    FastTree -gtr -nt $fasta > ${fasta.baseName}.nwk
    """
}
