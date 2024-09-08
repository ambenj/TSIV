#!/bin/bash
#SBATCH --job-name=MCP_fasta
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --account=kingsley
#SBATCH --partition=batch


# Inputs
LIST="/labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/06_blastn_iridoviruses/MCP_genes_curated.txt"
FASTA_DIR="/labs/kingsley/ambenj/TSIV/analysis/other_virus_sequences"
OUT_DIR="/labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/06_blastn_iridoviruses"

# Prepare output files
ALL_GENOMES="${OUT_DIR}/MCP_all_genomes.fasta"
CORE_TREE="${OUT_DIR}/MCP_coretree_genomes.fasta"
MEG="${OUT_DIR}/MCP_megalocytivirus_genomes.fasta"

echo $ALL_GENOMES
echo $CORE_TREE
echo $MEG

# Make directory if necessary and remove any prexisting files
mkdir -p $OUT_DIR
rm $ALL_GENOMES
rm $CORE_TREE
rm $MEG

module load samtools/1.20

# Remove carriage returns in input file
sed -i 's/\r//' $LIST

# Initialize a flag to skip the first line (header)
skip_header=true

# Loop through each line of the file
while IFS=$'\t' read -r virus_accession genus core seq
do
    # Skip the header line
    if [ "$skip_header" = true ]; then
        skip_header=false
        continue
    fi

    echo "Processing: $virus_accession $genus $core $seq"

    # Get relevant input fasta file
    FASTA="${FASTA_DIR}/${virus_accession}.1/${virus_accession}.1_orf_nucleotides.fasta"

    # Extract ORF sequence and put it in correct ORF file, renaming it after the virus_accession
    FASTA_ENTRY=$( samtools faidx $FASTA $seq | sed 's/^>'${seq}'/>'${virus_accession}'/' )
    echo "$FASTA_ENTRY" >> "$ALL_GENOMES"

    # Add entry to core genes tree fasta if it was used for the core gene analysis
    if [ "$core" == "yes" ]; then
        echo "$FASTA_ENTRY" >> "$CORE_TREE"
    fi

    # Add entry to megalocytivirus fasta if it is a megalocytivirus
    if [ "$genus" == "Megalocytivirus" ]; then
        echo "$FASTA_ENTRY" >> "$MEG"
    fi


done < "$LIST"

