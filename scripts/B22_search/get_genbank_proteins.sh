#!/bin/bash
#SBATCH --job-name=genbank
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --account=kingsley
#SBATCH --partition=batch


# Inputs
LIST=$1
OUT_FILE=$2

module load entrez-direct/11.0

# Remove carriage returns in input file
sed -i 's/\r//' $LIST

# Remove old outfile if it exists
rm $OUT_FILE

# Initialize a flag to skip the first line (header)
skip_header=true

# Loop through each line of the file
while IFS=$'\t' read -r accession title
do

    # Skip the header line
    if [ "$skip_header" = true ]; then
        skip_header=false
        continue
    fi

    echo "Processing: $accession $title"

    # Retrive protein sequences from genbank
    efetch -db protein -format fasta -id $accession >> $OUT_FILE
    
done < "$LIST"

