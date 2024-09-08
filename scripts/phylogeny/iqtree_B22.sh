#!/bin/bash
#SBATCH --job-name=iqtreeB
#SBATCH --ntasks=24
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Get input files
DIR="/labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/08_B22_phylogeny"
ALIGN_FILE="${DIR}/B22_MAFFT_alignment.fasta"
OUT_PREFIX="${DIR}/B22_MAFFT_alignment_boot100"

# Run IQ-TREE
## Default is using ModelFinderPlus to find best fit model
## Using standard bootstrap (minimum recommended is 100)
iqtree -s $ALIGN_FILE -b 100 -T AUTO -ntmax 24