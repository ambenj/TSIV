#!/usr/bin/env python
from os.path import join, abspath, expanduser
import re

# Read in config file
OUTDIR = config['out_dir']
FASTQ_DIR = config['fastq_dir']
F_EXT = config['f_ext']
R_EXT = config['r_ext']
SAMP = config['samp_prefixes']


# Convert out directory to absolute path
if OUTDIR[0] == '~':
    OUTDIR = expanduser(OUTDIR)
OUTDIR = abspath(OUTDIR)

###############################################################

rule all:
    input:
        expand(join(OUTDIR, "assembly/01_preprocess/pre_fastqc/{samp}_2_fastqc.html"), samp=SAMP),
        expand(join(OUTDIR, "assembly/02_spades/{samp}_M/quast/report.tsv"), samp=SAMP),
        expand(join(OUTDIR, "assembly/02_spades/{samp}_M/remapped/{samp}.aln.sorted_coverage.txt"), samp=SAMP),

###############################################################
################## Preprocess unmapped reads ##################
###############################################################

# quality check unmapped reads before processing
rule pre_fastqc:
    input:
        R1 = join(FASTQ_DIR, "{samp}" + F_EXT),
        R2 = join(FASTQ_DIR, "{samp}" + R_EXT),
    output:
        R1 = join(OUTDIR, "assembly/01_preprocess/pre_fastqc/{samp}_1_fastqc.html"),
        R2 = join(OUTDIR, "assembly/01_preprocess/pre_fastqc/{samp}_2_fastqc.html"),
    params:
        outdir = join(OUTDIR, "assembly/01_preprocess/pre_fastqc/")
    resources:
        mem = 50,
        time = 6
    shell: """
        module load fastqc/0.12.1

        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir}
    """

# Trim adaptor sequences and qc with fastqc
rule trim_paired:
    input:
        R1 = join(FASTQ_DIR, "{samp}" + F_EXT),
        R2 = join(FASTQ_DIR, "{samp}" + R_EXT),
    output:
        R1 = join(OUTDIR, "assembly/01_preprocess/{samp}_1_val_1.fq"),
        R2 = join(OUTDIR, "assembly/01_preprocess/{samp}_2_val_2.fq"),
        qc1 = join(OUTDIR, "assembly/01_preprocess/post_fastqc/{samp}_1_val_1_fastqc.html"),
        qc2 = join(OUTDIR, "assembly/01_preprocess/post_fastqc/{samp}_2_val_2_fastqc.html")
    params:
        outdir = join(OUTDIR, "assembly/01_preprocess/"),
        fastqc_dir = join(OUTDIR, "assembly/01_preprocess/post_fastqc/")
    resources:
        mem = 32,
        time = 12
    threads: 4
    shell: """
        module load trim_galore/0.6.10
        module load cutadapt/3.4
        mkdir -p {params.outdir}
        mkdir -p {params.fastqc_dir}
        trim_galore --quality 20 \
            --length 60 \
            --output_dir {params.outdir} \
            --paired {input.R1} {input.R2} \
            --fastqc_args "-t {threads} --outdir {params.fastqc_dir}"
    """

###############################################################
###################### Assemble reads  ########################
###############################################################

# Run SPAdes assembler using metagenomics pipeline
rule spades_assembly:
    input:
        R1 = rules.trim_paired.output.R1,
        R2 = rules.trim_paired.output.R2,
    output:
        join(OUTDIR, "assembly/02_spades/{samp}_M/contigs.fasta"),
    threads: 16
    params:
        outdir = join(OUTDIR, "assembly/02_spades/{samp}_M"),
    resources:
        mem = 48,
        time = 48,
    shell: """
        module load spades/3.15.2
        spades.py --meta -o {params.outdir} -m {resources.mem} -t {threads} -1 {input.R1} -2 {input.R2}
    """

# Evaluate assembly quality
rule spades_quast:
    input:
        rules.spades_assembly.output
    output:
        join(OUTDIR, "assembly/02_spades/{samp}_M/quast/report.tsv"),
    params:
        outdir = join(OUTDIR, "assembly/02_spades/{samp}_M/quast/"),
    resources:
        mem = 8,
        time = 1
    shell: """
        module load quast/5.0.2
        quast.py -o {params.outdir} {input} --fast
    """

###############################################################
############## Evaluate assembly coverage  ####################
###############################################################

# Map reads back to contigs to try to ID which contigs with consistent high coverage
rule map_to_contigs:
    input:
        contigs = rules.spades_assembly.output,
        R1 = rules.trim_paired.output.R1,
        R2 = rules.trim_paired.output.R2,
    output:
        bam = join(OUTDIR, "assembly/02_spades/{samp}_M/remapped/{samp}.aln.sorted.bam"),
    params:
        outdir = join(OUTDIR, "assembly/02_spades/{samp}_M/remapped/"),
    threads: 4
    resources:
        mem = 24,
        time = 48
    shell: """
        module load bwa/0.7.18
        module load samtools/1.20

        # Index contigs
        bwa index {input.contigs}

        # Map reads to contigs
        bwa mem -M -t {threads} {input.contigs} {input.R1} {input.R2} | \
        samtools sort -O bam -o {output.bam} -
        samtools index {output.bam}
    """

# Get coverage report after mapping reads back to contigs
rule contig_coverage:
    input:
        bam = rules.map_to_contigs.output.bam,
        contigs = rules.spades_assembly.output,
    output:
        wg = join(OUTDIR, "assembly/02_spades/{samp}_M/remapped/{samp}.aln.sorted_coverage.txt"),
    resources:
        time = 24,
        mem = 16
    shell: """
        module load samtools/1.20

        printf "contig\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n" > {output.wg}

        # Get coverage for whole genome
        samtools coverage -H -q 3 {input.bam} >> {output.wg}
    """