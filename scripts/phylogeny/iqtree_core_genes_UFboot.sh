#!/bin/bash
#SBATCH --job-name=iqtreeUF
#SBATCH --ntasks=24
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Get input files
DIR="/labs/kingsley/ambenj/TSIV/analysis/assembly/05_core_phylogeny"
ALIGN_FILE="${DIR}/24_core_proteins_alignment.phy"
OUT_PREFIX="${DIR}/24_core_proteins_alignment_UFboot1000"

# Run IQ-TREE
## Default is using ModelFinderPlus to find best fit model
## Using standard bootstrap (minimum recommended is 100)
iqtree -s $ALIGN_FILE -B 1000 -T AUTO -ntmax 24 --prefix $OUT_PREFIX