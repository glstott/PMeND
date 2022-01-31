#!/bin/bash
#SBATCH --job-name=nj_treegen
#SBATCH --partition=batch
#SBATCH --ntasks=32
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=nj_treegen.log
#SBATCH --mail-user=gs69042@uga.edu
#SBATCH --mail-type=END,FAIL

ml FastTree

input_dir=/scratch/gs69042/bahllab/data/

cd $input_dir
for file in `ls *aligned*fasta`; do
	FastTree -gtr -nt $file > ${file%.fa*}.nwk
done
