#!/bin/bash
#SBATCH --job-name=core
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --account=kingsley
#SBATCH --partition=batch


# Inputs
PROTEINS="/labs/kingsley/ambenj/TSIV/analysis/assembly/04_blastp_iridoviruses/core_genes_TSIV.txt"
PROTEIN_FASTA="/labs/kingsley/ambenj/TSIV/analysis/assembly/03_annotate/geneMarkS_TSIV_STMY_X_2011_03_M_contigs/TSIV_STMY_X_2011_03_M_contigs_GeneMarkS.gff.faa"
OUT_DIR="/labs/kingsley/ambenj/TSIV/analysis/assembly/04_blastp_iridoviruses/core_proteins"

# Make directory if necessary and remove any prexisting files
mkdir -p $OUT_DIR
rm ${OUT_DIR}/TSIV_core_*_proteins.fasta

module load samtools/1.20

# Remove carriage returns in input file
sed -i 's/\r//' $PROTEINS

# Initialize a flag to skip the first line (header)
skip_header=true

# Loop through each line of the file
while IFS=$'\t' read -r order gene_name orf
do
    # Skip the header line
    if [ "$skip_header" = true ]; then
        skip_header=false
        continue
    fi

    echo "Processing: $order $gene_name $orf"

    # Get relevant output files
    ORF_FASTA="${OUT_DIR}/TSIV_core_${order}_proteins.fasta"

    # Extract ORF sequence and put it in correct ORF file, renaming it after the virus_accession
    samtools faidx $PROTEIN_FASTA $orf | sed 's/^>'${orf}'/>TSIV/' >> ${ORF_FASTA}
    
done < "$PROTEINS"

