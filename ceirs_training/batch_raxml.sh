#!/bin/bash
#SBATCH --job-name=raml
#SBATCH --partition=batch
#SBATCH --ntasks=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00


ml raxml
raxmlHPC -s Asia2017H3N2_subset.mafftout.v3.fasta -n Asia2017H3N2_subset.mafftout.v3.phy -m GTRGAMMA -N autoMRE -p 456 -x 456 -T 8
