#!/bin/bash
#SBATCH --job-name=genbank
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH --account=kingsley
#SBATCH --partition=batch


# Inputs
LIST="/labs/kingsley/ambenj/TSIV/resources/240814_genbank_accession_list.txt"
FASTA_DIR="/labs/kingsley/ambenj/TSIV/analysis/other_virus_sequences"
SHORTCUT_DIR="/labs/kingsley/ambenj/TSIV/analysis/other_virus_sequences/protein_DBs"

cd $SHORTCUT_DIR

module load entrez-direct/11.0

# Remove carriage returns in input file
sed -i 's/\r//' $LIST


# Loop through each line of the file
while IFS=$'\t' read -r accession virus
do

    echo "Processing: $accession $virus"

    # Prepare directory
    FILE_PREFIX="${virus}_${accession}"
    mkdir -p ${FASTA_DIR}/${FILE_PREFIX}

    # Retrive protein sequences from genbank
    efetch -db nucleotide -format fasta_cds_aa -id $accession > ${FASTA_DIR}/${FILE_PREFIX}/${FILE_PREFIX}_orf_proteins.fasta

    # Make shortcut for blast search
    ln -s "../${FILE_PREFIX}/${FILE_PREFIX}_orf_proteins.fasta"
    
done < "$LIST"

