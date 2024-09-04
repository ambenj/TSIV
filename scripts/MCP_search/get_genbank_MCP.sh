#!/bin/bash
#SBATCH --job-name=genbank
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --account=kingsley
#SBATCH --partition=batch


# Inputs
OUT="/labs/kingsley/ambenj/TSIV/analysis/other_virus_sequences/taxid308906_MCP_nucleotides.fasta"

module load entrez-direct/11.0

# Retrive nucleotide sequences from genbank for Megalocytiviruses and major capsid protein
esearch -db nucleotide -query "txid308906[Organism] AND major capsid protein" | efetch -format fasta > $OUT

