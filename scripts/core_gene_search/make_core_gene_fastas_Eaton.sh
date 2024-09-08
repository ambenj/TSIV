#!/bin/bash
#SBATCH --job-name=core
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --account=kingsley
#SBATCH --partition=batch


# Inputs
PROTEINS="/labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/04_blastp_iridoviruses/core_genes_curated_long.txt"
FASTA_DIR="/labs/kingsley/ambenj/TSIV/analysis/other_virus_sequences"
OUT_DIR="/labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/04_blastp_iridoviruses/core_proteins"

# Make directory if necessary and remove any prexisting files
mkdir -p $OUT_DIR
rm ${OUT_DIR}/Eaton_core_*_proteins.fasta

module load samtools/1.20

# Remove carriage returns in input file
sed -i 's/\r//' $PROTEINS

# Initialize a flag to skip the first line (header)
skip_header=true

# Loop through each line of the file
while IFS=$'\t' read -r order gene_name virus_accession orf
do
    # Skip the header line
    if [ "$skip_header" = true ]; then
        skip_header=false
        continue
    fi

    echo "Processing: $order $gene_name $virus_accession $orf"

    # Get relevant input fasta file and output files
    PROTEIN_FASTA="${FASTA_DIR}/${virus_accession}.1/${virus_accession}.1_orf_proteins.fasta"
    ORF_FASTA="${OUT_DIR}/Eaton_core_${order}_proteins.fasta"

    # Extract ORF sequence and put it in correct ORF file, renaming it after the virus_accession
    samtools faidx $PROTEIN_FASTA $orf | sed 's/^>'${orf}'/>'${virus_accession}'/' >> ${ORF_FASTA}
    
done < "$PROTEINS"

