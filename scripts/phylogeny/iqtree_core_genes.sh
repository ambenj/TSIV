#!/bin/bash
#SBATCH --job-name=iqtreeB
#SBATCH --ntasks=16
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Get input files
DIR="/labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/05_core_phylogeny"
ALIGN_FILE="${DIR}/24_core_proteins_MAFFT_alignment.fasta"
OUT_PREFIX="${DIR}/24_core_proteins_MAFFT_alignment_boot100"

# Run IQ-TREE
## Default is using ModelFinderPlus to find best fit model
## Using standard bootstrap (minimum recommended is 100)
iqtree -s $ALIGN_FILE -b 100 -T AUTO -ntmax 16 --prefix $OUT_PREFIX