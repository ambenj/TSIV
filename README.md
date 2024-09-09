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

## Assembly and annotation of TSIV genome
### Assembly
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

### Annotation of the TSIV genome
1. Calculate GC content
2. Predict ORFs using GeneMarkS
3. BLAST ORF nucleotide sequences against the nt database, ISKNV, and FV3
4. BLAST ORF protein sequences against the nr database, ISKNV, and FV3

```bash
module load snakemake
snakemake --configfile configs/240905_annotate_TSIV_STMY_2012_42.yaml --snakefile scripts/assembly_and_annotation/annotate.smk --profile scg --jobs 300 --restart-times 0 --rerun-incomplete
```
