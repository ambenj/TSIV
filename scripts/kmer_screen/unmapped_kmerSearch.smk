#!/usr/bin/env python
from os.path import join, abspath, expanduser
import re

# Read in config file
OUTDIR = config['out_dir']
BAM_DIR = config['bam_dir']
EXT = config['ext']
SAMP = config['samp_prefixes']
KMER_LENGTH = config['kmer_length'] 
KMER_DIFFS = config['kmer_diffs']

# Convert out directory to absolute path
if OUTDIR[0] == '~':
    OUTDIR = expanduser(OUTDIR)
OUTDIR = abspath(OUTDIR)

# TSIV sequences
ATPase="CCGTGCTAATCAAATCTATCATTGCGGCCAAGCGCCACATCATACCGGCCGCCGTGGTCATCTCGGGCTCCGAGGAGGCCAACCACTTCTACAGCAAGCTGCTGCCCCAATGCTTTGTGTACTCAAAGTTTGACCCCGACATCATCATGCGCGTTAAGCAGCGCCAGCTGGCCCTGAAAAACATCGACCCTGAGCACTCGTGGCTGCTGCTCGTGTTTGACGACTGCATGGACAACGCCAAGGTCTTCAACCACGAGTCTGTCATGGACCTGTTCAAGAACGGGCGCCACTGGAACGTGCTGGTCATCATCGCGTCGCAGTACATAATGGACCTGAACGCGAGTCTGAGGTGCTGCATAGACGGCATATTTCTATTCACTGAGACAAGCCAGACATGTGTGGACAAGATATACAAACAGTTTGGTGGCAACATACCCAAACAGACATTTCACACCCTCATGGACAAGGTCACGCAGGACCACACGTGCCTCTACATAGACAACACCACGTCCCGGCAAAGGTGGGAGGACATGGTGCGCTACTACAAGGCGCCGCTCCTGGCCGACGGCGACGTGGCATTTGGCTTTGCCGACTACAAGGCGCACGGCGCCTCTACGATGTAAAGAGGTGCACGGCGCCCTTCTGCAGCGAGTGCTCGGGGTCGATGCTGCGCGTCAGAAACAGACCCTCTGTGCGCGGCTTGATGGGCACATCCCCG"
MCP="GGCTTCCAGAAAAAAACAAAAAACTAATATATAAAACAACAGAATGTCCTCAATCTCGGGTGCCAATGTGACGAGCGGGTTCATCGACATCGCCGCCTACGATGCCATGGAGACGCATCTGTACGGCGGCGACCACGCCGTGACATACTTTACCCGCGAGACCGTGCGCAGCTCCTGGTACAGCAAGCTGCCCGTCACACTGTCCAAGCAGACCGGGCACGCCGACTTTGGGCACGAATTCAGCGTGACCGTGGCGCGAGGCGGCGACTACCTGCTCAACGTGTGGCTGCGTGTGAAGCTGCCCGTCATCACCGCCAGCAAGGAGAACGCCCGCGTGCGCTGGTGCGACAACTTTATGCACAATCTGGTGGAGGAGGTGTCTCTGTCCTTCAACGACCTGGTGGCGCAGACCATCACCAGCGAGTTCCTGGACTTCTGGAGCACCTGCATGATGCCTGGCAGCAAGCAGTCCGGCTACAACAAGATGATTGGCATGCGCAGCGACCTGGTGACAGGCGCCACTACTGGCCAGATCCTGCCCGGCCACTACCTCAACCTGCCCATCCCGCTGTTCTTCACCAGGGACACAGGCCTGGCGCTGCCCACCGTGTCTCTGCCCTACAACGAGGTGCGCATCCACTTCAAGTTGAGGAAGTGGGAGGATCTGCTGATTGCACAGTCCACCCAGGACGCGTTCACCATCGAGACGGCCCAACTTGCTGATATCACCAACGTCACCCCGTCGCTCAGCAACGTCTCTGTGATGGGCACATACGCCATCCTGACCAGCGAGGAGCGCGAGACCGTGGCCATGTCCAGCCGCACCATGCTCATCGAGCAGTGTCAGGTGGCACCCCGCGTGCCCGTGTCCGCCGCTGACACCTCACTGGTGCACCTGGACTTGCGTTTCAGCCACCCCGTCAAGGCCCTGTTCTTTGCCGTCAAGAATGTGACACATAATAACGTGCATAGCAACTACACAGCGGCCAGCCCCGTGTACGTAAACCCCAACGTCATCCTGCCGGCACACGCCACCAACCCGCTGTCGGAGGTGTCGCTCATCTACGAGAACACCGCCAGGCTGCACCAGATGGGCGTCGACTACTTCACCTCGGTGGACCCCTACTACTTTGCGCCCAGTATGCCCGAGATGGACGGTGTCATGGTGTACTGCTACACCCTGGACATGGGCAGCATCAACCCCATGGGCTCCACCAACTACGGCCGCCTGTCCAACGTGTCGCTGGCCTGCAAAGTGTCAGAGAACGCCCGCACCACAGCCGCTGGAGGCGGTGCAAACGGCAGCGGCTACACCGTGCCGCAGAAGTTTGAACTGGTGGTGATCGCTGTCAACCACAACATTATGAAGATTGCCAACGGCGCCG"
DNAPol="CAGCATCATCATCTCAAAAAACATATGCTACACCACGCTGGTGGACAGCGGCGGCGAGGAGTATGCCTGGCAGGAGCACCAGGGCTGTGAGCACGACCCAGCGCACGCGCAGCGCGAGGCCCTGGGCAGTGAGATAGGCGCGCTGCAGTGTGCCATTGCGGCGCTACCGCGCAAGGCCACGCAGGAGCGTGCGCGCCTCCGGGAACGCGTTGCCGAGATGAAGGTCCGCCACTCGAGCATGGCGCCGGCCTCCGTCAAGTGCAACGTGTTCAACTACAGGTTCACGCGCGAGCGCGAGGGCGTGCTGCCGCGTGTGCTGCGCAACCTGCTGGACAGCAGGGCGGACATACGCGAGCGCATGAAGGGCCTCTCGGACGTAGACATGCTCGGCGTCATGGACAAGAGGCAGCTGGCGTACAAAATCAGCGCCAACTCGGTGTACGGCACCATGGGCACACAGCGTGGCTACCTGCCCTTCATGGCGGGCGCCATGACCACCACCTACTGCGGACGCAAGCTCATCGAAAAGGCGGCGGAGCTGCTGCGCACCGTAGTGGGCGCAACCCTGGTG"


###############################################################

rule all:
    input:
        join(OUTDIR, "01_unmapped_reads/fastqc/multiqc_report.html"),
        expand(join(OUTDIR, "02_kmer_search/DNAPol/{samp}_stats.txt"), samp=SAMP),

###############################################################
######### Extract and process unmapped reads ##################
###############################################################

# Extract unmapped reads and convert to fastq
rule extract:
    input:
        bam = join(BAM_DIR, "{samp}" + EXT)
    output:
        R1 = join(OUTDIR, "01_unmapped_reads/{samp}_1.fq"),
        R2 = join(OUTDIR, "01_unmapped_reads/{samp}_2.fq"),
        s = join(OUTDIR, "01_unmapped_reads/{samp}_s.fq"),
    threads: 8
    resources:
        mem = 24,
        time = 24
    shell: """
        module load samtools/1.13
        cd $TMPDIR
        samtools sort -n -@ {threads} {input.bam} | \
        samtools fastq -@ {threads} -t -T BX -f 4 -1 {output.R1} -2 {output.R2} -s {output.s} -
    """

# quality check unmapped reads
rule fastqc:
    input:
        R1 = join(OUTDIR, "01_unmapped_reads/{samp}_1.fq"),
        R2 = join(OUTDIR, "01_unmapped_reads/{samp}_2.fq"),
        s = join(OUTDIR, "01_unmapped_reads/{samp}_s.fq"),
    output:
        R1 = join(OUTDIR, "01_unmapped_reads/fastqc/{samp}_1_fastqc.html"),
        R2 = join(OUTDIR, "01_unmapped_reads/fastqc/{samp}_2_fastqc.html"),
        s = join(OUTDIR, "01_unmapped_reads/fastqc/{samp}_s_fastqc.html"),
    params:
        outdir = join(OUTDIR, "01_unmapped_reads/fastqc/")
    resources:
        mem = 48,
        time = 12
    shell: """
        module load fastqc/0.12.1

        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir}
    """

# compile fastqc files into single report
rule multiqc:
    input:
        expand(join(OUTDIR, "01_unmapped_reads/fastqc/{samp}_{read}_fastqc.html"), samp=SAMP, read=["1","2","s"])
    output:
        join(OUTDIR, "01_unmapped_reads/fastqc/multiqc_report.html")
    params:
        dir = join(OUTDIR, "01_unmapped_reads/fastqc/")
    shell:"""
        module load multiqc/1.9
        multiqc --force {params.dir} -o {params.dir}
    """

###############################################################
######################## Search kmers #########################
###############################################################

rule search_kmers:
    input:
        R1 = join(OUTDIR, "01_unmapped_reads/{samp}_1.fq"),
        R2 = join(OUTDIR, "01_unmapped_reads/{samp}_2.fq"),
    output:
        R1_ATPase = join(OUTDIR, "02_kmer_search/ATPase/{samp}_1.fq"),
        R2_ATPase = join(OUTDIR, "02_kmer_search/ATPase/{samp}_2.fq"),
        stats_ATPase = join(OUTDIR, "02_kmer_search/ATPase/{samp}_stats.txt"),
        R1_MCP = join(OUTDIR, "02_kmer_search/MCP/{samp}_1.fq"),
        R2_MCP = join(OUTDIR, "02_kmer_search/MCP/{samp}_2.fq"),
        stats_MCP = join(OUTDIR, "02_kmer_search/MCP/{samp}_stats.txt"),
        R1_DNAPol = join(OUTDIR, "02_kmer_search/DNAPol/{samp}_1.fq"),
        R2_DNAPol = join(OUTDIR, "02_kmer_search/DNAPol/{samp}_2.fq"),
        stats_DNAPol = join(OUTDIR, "02_kmer_search/DNAPol/{samp}_stats.txt"),
    shell:"""
        module load bbmap/39.01

        # Search for TSIV ATPase
        bbduk.sh in={input.R1} in2={input.R2} outm={output.R1_ATPase} outm2={output.R2_ATPase} literal={ATPase} k={KMER_LENGTH} hdist={KMER_DIFFS} stats={output.stats_ATPase}

        # Search for TSIV MCP
        bbduk.sh in={input.R1} in2={input.R2} outm={output.R1_MCP} outm2={output.R2_MCP} literal={MCP} k={KMER_LENGTH} hdist={KMER_DIFFS} stats={output.stats_MCP}

        # Search for TSIV DNAPol
        bbduk.sh in={input.R1} in2={input.R2} outm={output.R1_DNAPol} outm2={output.R2_DNAPol} literal={DNAPol} k={KMER_LENGTH} hdist={KMER_DIFFS} stats={output.stats_DNAPol}

    """