#!/usr/bin/env python
from os.path import join, abspath, expanduser

# Read in config file
OUTDIR = config['out_dir']
QUERY_FA = config['query_nucl_fa']
SUBJECT_FA_DIR = config['subject_fa_dir']
SUBJECT_FA_PREFIXES = config['subject_nucl_fa_prefixes']


# Convert out directory to absolute path
if OUTDIR[0] == '~':
    OUTDIR = expanduser(OUTDIR)
OUTDIR = abspath(OUTDIR)

###############################################################

rule all:
    input:
        expand(join(OUTDIR, "{subject_prefix}_blastn.tsv"), subject_prefix=SUBJECT_FA_PREFIXES),

###############################################################

# Generate a blastn database for input CDS sequences
rule create_blastn_database:
    input:
        nucleotide_fa = join(SUBJECT_FA_DIR, "{subject_prefix}.fasta"),
    output:
        nucleotide_db = join(SUBJECT_FA_DIR, "{subject_prefix}.fasta.ndb"),
    resources:
        mem = 8,
        time = 1
    shell: """
        module load ncbi-blast/2.15.0+
        makeblastdb -in {input} -dbtype nucl
    """


# Blast input query sequences against each custom database
rule blastn:
    input:
        query = {QUERY_FA},
        subject_db = rules.create_blastn_database.output.nucleotide_db,
        nucleotide_fa = join(SUBJECT_FA_DIR, "{subject_prefix}.fasta"),
    output:
        blast = join(OUTDIR, "{subject_prefix}_blastn.tsv"),
    threads: 4
    resources:
        mem = 16,
        time = 4
    shell: """
        module load ncbi-blast/2.15.0+
        blastn -task dc-megablast -query {input.query} -db {input.nucleotide_fa} -num_threads {threads} -out {output.blast} -evalue 1 \
        -outfmt "6 qseqid sseqid evalue stitle pident length mismatch gapopen qstart qend sstart send bitscore"
    """