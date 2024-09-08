#!/bin/bash
#SBATCH --job-name=sdt_muscle
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --mem=16G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Run sequence demarcation tool with mafft setting for pairwise identity calculations

# Note: need to run this in environment with python 2.7
# conda activate py27
# sbatch SDT_mafft.sh <input_fasta> <out_dir>

FASTA=$1  # Absolute path to input fasta file
OUTDIR=$2  # Path to outdir
SDT="/labs/kingsley/ambenj/tools/SDT_Linux64"

# Make outdir if it does not exist and navigate to it
mkdir -p ${OUTDIR}
cd $OUTDIR
mkdir -p output

# Run SDT (updated with paths on system)
python ${SDT}/SDT_Linux64_paths.py $FASTA muscle