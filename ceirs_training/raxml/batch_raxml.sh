#!/bin/bash
#SBATCH --job-name=raxml
#SBATCH --partition=batch
#SBATCH --ntasks=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00


ml raxml
raxmlHPC -s /home/gs69042/CEIRS-Training-Taiwan-2019/CIERStraining-Asia2017HA/RAxML/Asia2017H3N2_subset.mafftout.final.fasta -n Asia2017H3N2_subset.mafftout.v3.phy -m GTRGAMMA -p 456 -N 3

