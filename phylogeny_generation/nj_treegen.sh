#!/bin/bash
#SBATCH --job-name=mafft_fasttree
#SBATCH --partition=batch
#SBATCH --ntasks=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=mafft_fasttree.log
#SBATCH --mail-user=gs69042@uga.edu
#SBATCH --mail-type=END,FAIL

ml MAFFT
ml FastTree

input_dir=/scratch/gs69042/bahllab/data/

mafft --auto $input_dir/GISAID_combined.fasta > $input_dir/GISAID_combined.aligned.fasta

