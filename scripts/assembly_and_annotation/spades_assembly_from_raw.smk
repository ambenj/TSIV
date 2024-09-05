#!/usr/bin/env python
from os.path import join, abspath, expanduser
import re

# Read in config file
OUTDIR = config['out_dir']
REF=config['ref']
FASTQ_DIR = config['fastq_dir']
F_EXT = config['f_ext']
R_EXT = config['r_ext']
SAMP = config['samp_prefixes']
SUBSAMP = config['subsamp']

# Convert out directory to absolute path
if OUTDIR[0] == '~':
    OUTDIR = expanduser(OUTDIR)
OUTDIR = abspath(OUTDIR)

###############################################################

rule all:
    input:
        expand(join(OUTDIR, "assembly/01_preprocess/pre_fastqc/{samp}_R2_001_fastqc.html"), samp=SAMP),
        expand(join(OUTDIR, "assembly/01_preprocess/{samp}.sorted_coverage.txt"), samp=SAMP),
        expand(join(OUTDIR, "assembly/02_spades/{samp}_down{subsamp}/quast/report.tsv"), samp=SAMP, subsamp=SUBSAMP),
        expand(join(OUTDIR, "assembly/02_spades/{samp}_down{subsamp}/remapped/{samp}.aln.sorted_coverage.txt"), samp=SAMP, subsamp=SUBSAMP),

###############################################################
##################### Preprocess raw reads ####################
###############################################################

# quality check unmapped reads before processing
rule pre_fastqc:
    input:
        R1 = join(FASTQ_DIR, "{samp}" + F_EXT),
        R2 = join(FASTQ_DIR, "{samp}" + R_EXT),
    output:
        R1 = join(OUTDIR, "assembly/01_preprocess/pre_fastqc/{samp}_R1_001_fastqc.html"),
        R2 = join(OUTDIR, "assembly/01_preprocess/pre_fastqc/{samp}_R2_001_fastqc.html"),
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
        R1 = join(OUTDIR, "assembly/01_preprocess/{samp}_R1_001_val_1.fq.gz"),
        R2 = join(OUTDIR, "assembly/01_preprocess/{samp}_R2_001_val_2.fq.gz"),
        qc1 = join(OUTDIR, "assembly/01_preprocess/post_fastqc/{samp}_R1_001_val_1_fastqc.html"),
        qc2 = join(OUTDIR, "assembly/01_preprocess/post_fastqc/{samp}_R2_001_val_2_fastqc.html")
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
            --length 200 \
            --output_dir {params.outdir} \
            --paired {input.R1} {input.R2} \
            --fastqc_args "-t {threads} --outdir {params.fastqc_dir}"
    """

###############################################################
########## Remove stickleback reads downsample ################
###############################################################

rule remove_stickleback:
    input:
        R1 = rules.trim_paired.output.R1,
        R2 = rules.trim_paired.output.R2,
    output:
        bam = join(OUTDIR, "assembly/01_preprocess/{samp}.sorted.bam"),
        R1 = join(OUTDIR, "assembly/01_preprocess/{samp}_ref_unmapped_1.fq"),
        R2 = join(OUTDIR, "assembly/01_preprocess/{samp}_ref_unmapped_2.fq"),
        s = join(OUTDIR, "assembly/01_preprocess/{samp}_ref_unmapped_s.fq"),
    threads: 8
    resources:
        mem = 24,
        time = 24
    shell:"""
        module load bwa/0.7.18
        module load samtools/1.20

        # Index contigs
        bwa index {REF}

        # Map reads to contigs
        bwa mem -M -t {threads} {REF} {input.R1} {input.R2} | \
        samtools sort -O bam -o {output.bam} -
        samtools index {output.bam}

        # Extract unmapped reads
        samtools sort -n -@ {threads} {output.bam} | \
        samtools fastq -@ {threads} -t -T BX -f 4 -1 {output.R1} -2 {output.R2} -s {output.s} -
    """

# Get coverage of reads mapped to host genome
rule ref_coverage:
    input:
        bam = rules.remove_stickleback.output.bam,
    output:
        wg = join(OUTDIR, "assembly/01_preprocess/{samp}.sorted_coverage.txt"),
    resources:
        time = 24,
        mem = 8
    shell: """
        module load samtools/1.20

        printf "contig\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n" > {output.wg}

        # Get coverage for whole genome
        samtools coverage -H -q 3 {input.bam} >> {output.wg}
    """


# Dowsample reads for SPADES assembly
rule downsample_unmapped:
    input:
        R1 = rules.remove_stickleback.output.R1,
        R2 = rules.remove_stickleback.output.R2,
    output:
        R1 = join(OUTDIR, "assembly/01_preprocess/{samp}_ref_unmapped_down{subsamp}_1.fq"),
        R2 = join(OUTDIR, "assembly/01_preprocess/{samp}_ref_unmapped_down{subsamp}_2.fq"),
    resources:
        mem = 12,
        time = 12,
    shell: """
        module load seqtk/1.4-r130

        # Downsample reads
        seqtk sample -s100 {input.R1} {wildcards.subsamp} > {output.R1}
        seqtk sample -s100 {input.R2} {wildcards.subsamp} > {output.R2}

    """

###############################################################
###################### Assemble reads  ########################
###############################################################

# Run SPAdes assembler using metagenomics pipeline
rule spades_assembly:
    input:
        R1 = rules.downsample_unmapped.output.R1,
        R2 = rules.downsample_unmapped.output.R2,
    output:
        join(OUTDIR, "assembly/02_spades/{samp}_down{subsamp}/contigs.fasta"),
    threads: 8
    params:
        outdir = join(OUTDIR, "assembly/02_spades/{samp}_down{subsamp}"),
    resources:
        mem = 400,
        time = 72,
    shell: """
        module load spades/3.15.2

        # Run spades
        #spades.py --meta -o {params.outdir} -m {resources.mem} -t {threads} -1 {input.R1} -2 {input.R2}
        spades.py -o {params.outdir} -m {resources.mem} -t {threads} -1 {input.R1} -2 {input.R2}


        # Restart spades after memory problem
        # spades.py --restart-from last -o {params.outdir} -m {resources.mem} -t {threads}

    """

# Evaluate assembly quality
rule spades_quast:
    input:
        rules.spades_assembly.output
    output:
        join(OUTDIR, "assembly/02_spades/{samp}_down{subsamp}/quast/report.tsv"),
    params:
        outdir = join(OUTDIR, "assembly/02_spades/{samp}_down{subsamp}/quast/"),
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
        R1 = rules.remove_stickleback.output.R1,
        R2 = rules.remove_stickleback.output.R2,
    output:
        bam = join(OUTDIR, "assembly/02_spades/{samp}_down{subsamp}/remapped/{samp}.aln.sorted.bam"),
    params:
        outdir = join(OUTDIR, "assembly/02_spades/{samp}_down{subsamp}/remapped/"),
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
        wg = join(OUTDIR, "assembly/02_spades/{samp}_down{subsamp}/remapped/{samp}.aln.sorted_coverage.txt"),
    resources:
        time = 24,
        mem = 16
    shell: """
        module load samtools/1.20

        printf "contig\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n" > {output.wg}

        # Get coverage for whole genome
        samtools coverage -H -q 3 {input.bam} >> {output.wg}
    """