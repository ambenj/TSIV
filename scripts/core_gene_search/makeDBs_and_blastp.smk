#!/usr/bin/env python
from os.path import join, abspath, expanduser

# Read in config file
OUTDIR = config['out_dir']
QUERY_FA = config['query_prot_fa']
SUBJECT_FA_DIR = config['subject_fa_dir']
SUBJECT_FA_PREFIXES = config['subject_prot_fa_prefixes']


# Convert out directory to absolute path
if OUTDIR[0] == '~':
    OUTDIR = expanduser(OUTDIR)
OUTDIR = abspath(OUTDIR)

###############################################################

rule all:
    input:
        expand(join(OUTDIR, "{subject_prefix}_blastp.tsv"), subject_prefix=SUBJECT_FA_PREFIXES),

###############################################################

# Generate a BLASTp database for input protein sequences
rule create_blastp_database:
    input:
        protein_fa = join(SUBJECT_FA_DIR, "{subject_prefix}.fasta"),
    output:
        protein_db = join(SUBJECT_FA_DIR, "{subject_prefix}.fasta.pdb"),
    resources:
        mem = 8,
        time = 1
    shell: """
        module load ncbi-blast/2.15.0+
        makeblastdb -in {input} -dbtype prot
    """


# Blast input query sequences against each custom database
rule blastp:
    input:
        query = {QUERY_FA},
        subject_db = rules.create_blastp_database.output.protein_db,
        protein_fa = join(SUBJECT_FA_DIR, "{subject_prefix}.fasta"),
    output:
        blast = join(OUTDIR, "{subject_prefix}_blastp.tsv"),
    threads: 4
    resources:
        mem = 16,
        time = 4
    shell: """
        module load ncbi-blast/2.15.0+
        blastp -query {input.query} -db {input.protein_fa} -num_threads {threads} -out {output.blast} -evalue 1 \
        -outfmt "6 qseqid sseqid evalue stitle pident length mismatch gapopen qstart qend sstart send bitscore"
    """