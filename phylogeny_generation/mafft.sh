#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH --partition=batch
#SBATCH --ntasks=32
#SBATCH --mem=120G
#SBATCH --time=72:00:00
#SBATCH --output=mafft.log

ml MAFFT

input_dir=/scratch/gs69042/bahllab/data/

cd $input_dir
for file in `ls *-*.fasta`; do
	mafft --6merpair --thread 32 --addfragments $file reference.fasta > ${file%.fa*}.aligned.fasta
done
