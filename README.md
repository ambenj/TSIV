# TSIV Genome Analysis

This repository contains the code used for following paper:

Yoxsimer, A.M.; Offenberg, E.G.; Katzer, A.W.; Bell, M.A.; Massengill, R.L.; Kingsley, D.M. Genomic sequence of the threespine stickleback iridovirus (TSIV) from wild Gasterosteus in Stormy Lake, Alaska. 2024.

*Note:these scripts were configured for analysis on the SCG Cluster at Stanford University and most would require adjusting paths and software in order to run elsewhere.*

## Screen for kmers in stickleback genomes
1. Extract unmapped reads
2. Perform kmer search using bbduk.sh

```bash
# Run kmer search
module load snakemake 
snakemake --configfile configs/240724_unmapped_kmerSearch_206genomes.yaml --snakefile scripts/kmer_screen/unmapped_kmerSearch.smk --profile scg --jobs 300 --restart-times 0 --rerun-incomplete
```

## Assembly of the TSIV genome
1. Quality check raw reads
2. Trim reads
3. Map reads to stickleback reference genome
4. Downsample reads
4. Extract unmapped reads
5. Assemble reads using SPADES
6. Evaluate SPADES assembly with coverage

```bash
# Run SPADES assembly pipeline
module load snakemake
snakemake --configfile configs/240904_spades_assembly_from_raw_STMY_2012_42.yaml --snakefile scripts/assembly_and_annotation/spades_assembly_from_raw.smk --profile scg --jobs 300 --restart-times 0 --rerun-incomplete
```

```bash
# Extract likely TSIV contigs from best assembly
module load samtools/1.20
faidx contigs.fasta NODE_1_length_79043_cov_150.476925 > TSIV_STMY_2012_42.fasta
faidx contigs.fasta NODE_2_length_36080_cov_145.724195 >> TSIV_STMY_2012_42.fasta
```
* `NODE_1_length_79043_cov_150.476925` was renamed to `contig1`
* `NODE_2_length_36080_cov_145.724195` was renamed to `contig2`

## Annotation of the TSIV genome

### ORF prediction
1. Calculate GC content
2. Predict ORFs using GeneMarkS
3. BLAST ORF nucleotide sequences against the nt database, ISKNV, and FV3
4. BLAST ORF protein sequences against the nr database, ISKNV, and FV3

```bash
# Run annotation pipeline
module load snakemake
snakemake --configfile configs/240905_annotate_TSIV_STMY_2012_42.yaml --snakefile scripts/assembly_and_annotation/annotate.smk --profile scg --jobs 300 --restart-times 0 --rerun-incomplete
```

### Iridovirus core gene identification
1. BLAST TSIV ORF proteins against iridovirus protein databases
2. Curate annotations

```bash
# Run BLAST database pipeline
snakemake --configfile configs/240905_makeDBs_and_blastp.yaml --snakefile scripts/core_gene_search/makeDBs_and_blastp.smk --profile scg --jobs 300 --restart-times 0 --rerun-incomplete
```

Ran the following two scripts in R to analysis BLASTP hits and curate annotation tables:
* `scripts/Rscripts/240905_TSIV_contigs_annotation_analysis.Rmd`
* `scripts/Rscripts/240905_TSIV_STMY_2012_42_homolog_analysis.Rmd`

## Phylogenetic analysis and sequence identity matrices
### Core genes
1. Gather core genes into fasta files for alignment
2. Generate alignments in Geneious
3. Run phylogenetic analysis
4. Calculate pairwise sequence identity matrix

```bash
# Gather core genes into fasta files
sbatch scripts/core_gene_search/make_core_gene_fastas_Eaton.sh
sbatch scripts/core_gene_search/make_core_gene_fastas_Genbank.sh
sbatch scripts/core_gene_search/make_core_gene_fastas_TSIV.sh

# Join all protein sequences into one fasta file for each core gene
for n in {1..26}; do cat *_${n}_proteins.fasta > all_core_${n}_proteins.fasta; done
```

```bash
# After generating alignments in Geneious, run IQ-TREE to perform maximum likelihood analysis
conda activate tree
sbatch scripts/phylogeny/iqtree_core_genes.sh
```

```bash
# Prepare sequence identity matrx
conda activate py27
sbatch scripts/pairwise_identity/SDT_mafft.sh /labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/09_pairwise_identity/24_core_protein/24_core_proteins_concat.fasta /labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/09_pairwise_identity/24_core_protein/mafft
```

### Major capsid protein
1. BLAST TSIV ORF nucleotide sequences against iridovirus nucleotide databases
2. Curate list of MCP sequences
2. Curate MCP list in Geneious
3. Calculate pairwise sequence identity matrix

```bash
# Run nucleotide blast search pipeline
module load snakemake
snakemake --configfile configs/240905_makeDBs_and_blastn.yaml --snakefile scripts/MCP_search/makeDBs_and_blastn.smk --profile scg --jobs 300 --restart-times 0 --rerun-incomplete
```

Ran the following script in R to analysis nucleotide BLAST hits and curate annotation tables:
* `scripts/Rscripts/240905_TSIV_MCP_nucleotide_analysis.Rmd`

```bash
# Gather nucleotide sequences into fasta file
sbatch scripts/MCP_search/make_MCP_fasta.sh
```

```bash
# After curating sequences in Geneious, prepare sequence identity matrx
conda activate py27
sbatch scripts/pairwise_identity/SDT_mafft.sh /labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/09_pairwise_identity/MCP_nucleotide/MCP_nucleotides.fasta /labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/09_pairwise_identity/MCP_nucleotide/mafft
```

### B22 phylogeny
1. Curate B22 list
2. Generate alignments in Geneious
3. Run phylogenetic analysis

```bash
# After generating alignments in Geneious, run IQ-TREE to perform maximum likelihood analysis
conda activate tree
sbatch scripts/phylogeny/iqtree_B22.sh
```

## Collinearity Analysis
1. Determine arrangement of TSIV genome contigs by finding read pairs spanning contig1 and contig2
2. Curate core gene annotations and rearrange genomes for consistent collinear analysis plots
3. Prepare collinearity plots

```bash
# Extract discordant read pairs
module load samtools/1.20
samtools view -f 1 -F 12 STMY12-42.aln.sorted.bam | awk '$3 == "NODE_1_length_79043_cov_150.476925" && $7 == "NODE_2_length_36080_cov_145.724195" || $3 == "NODE_2_length_36080_cov_145.724195" && $7 == "NODE_1_length_79043_cov_150.476925"' > discordant_pairs.txt

# Plot read positions
Rscript scripts/Rscripts/240907_discordant_pairs_figure.R
```

Genome arrangements and core gene annotations were curated in the following script: `scripts/Rscripts/240906_prep_annotations_for_collinearity.Rmd`

```bash
# Make collinearity plot
Rscript scripts/Rscripts/240906_collinearity_figure.R
```