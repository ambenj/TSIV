#!/usr/bin/env python
from os.path import join, abspath, expanduser

# Read in config file
OUTDIR = config['out_dir']
ASSEMBLY = config['assembly']
SAMP = config['prefix']


# Convert out directory to absolute path
if OUTDIR[0] == '~':
    OUTDIR = expanduser(OUTDIR)
OUTDIR = abspath(OUTDIR)

###############################################################

rule all:
    input:
        expand(join(OUTDIR, "geneMarkS_{samp}/{samp}_GC.tsv"), samp=SAMP),
        expand(join(OUTDIR, "geneMarkS_{samp}", "{samp}_GeneMarkS_blastn.tsv"), samp=SAMP),
        expand(join(OUTDIR, "geneMarkS_{samp}", "{samp}_GeneMarkS_blastp.tsv"), samp=SAMP),
        expand(join(OUTDIR, "geneMarkS_{samp}", "{samp}_GeneMarkS_blastn_FV3_AY548484.1.tsv"), samp=SAMP),
        expand(join(OUTDIR, "geneMarkS_{samp}", "{samp}_GeneMarkS_blastp_FV3_AY548484.1.tsv"), samp=SAMP),


###############################################################
######### GeneMarkS for coding sequence prediction ############
###############################################################

# Run geneMarkS to annotate predicted coding sequences
rule geneMarkS:
    input:
        contigs = {ASSEMBLY}
    output:
        gff = join(OUTDIR, "geneMarkS_{samp}", "{samp}_GeneMarkS.gff"),
        faa = join(OUTDIR, "geneMarkS_{samp}", "{samp}_GeneMarkS.gff.faa"),
        fnn = join(OUTDIR, "geneMarkS_{samp}", "{samp}_GeneMarkS.gff.fnn"),
    resources:
        mem = 16,
        time = 2
    params:
        dir = join(OUTDIR, "geneMarkS_{samp}"),
    shell: """
        module load genemarks/4.3
        mkdir -p {params.dir}
        cd {params.dir}
        gmsn.pl --virus --format GFF -fnn -faa --output {wildcards.samp}_GeneMarkS.gff {input.contigs}
    """

# Blast CDS nucleotide sequences called by geneMarkS
rule blastn_genemark:
    input:
        rules.geneMarkS.output.fnn
    output:
        blast = join(OUTDIR, "geneMarkS_{samp}/{samp}_GeneMarkS_blastn.tsv"),
    params:
        db = "/reference/ncbi/blast/db/v5/2024-03-08/nt",
    threads: 16
    resources:
        mem = 48,
        time = 24
    shell: """
        module load ncbi-blast/2.15.0+
        blastn -task dc-megablast \
        -outfmt "6 qseqid sseqid evalue stitle pident length mismatch gapopen qstart qend sstart send bitscore sscinames scomnames staxid sacc" \
        -query {input} -db {params.db} -num_threads {threads} -out {output.blast} -evalue 1
    """

###############################################################
#################### Calculate CG content #####################
###############################################################

rule GC:
    input: {ASSEMBLY}
    output: join(OUTDIR, "geneMarkS_{samp}/{samp}_GC.tsv")
    shell:"""
        module load seqkit/2.3.1
        seqkit stats -a {input} > {output}
    """


###############################################################
######################## BLAST searches #######################
###############################################################

# Blast translated CDS called by geneMarkS
rule blastp_genemark:
    input:
        rules.geneMarkS.output.faa
    output:
        blast = join(OUTDIR, "geneMarkS_{samp}/{samp}_GeneMarkS_blastp.tsv"),
    params:
        db = "/reference/ncbi/blast/db/v5/2024-03-08/nr",
    threads: 32
    resources:
        mem = 48,
        time = 72
    shell: """
        module load ncbi-blast/2.15.0+
        blastp -query {input} -db {params.db} -num_threads {threads} -out {output.blast} -evalue 1 \
        -outfmt "6 qseqid sseqid evalue stitle pident length mismatch gapopen qstart qend sstart send bitscore sscinames scomnames staxid sacc"
    """

# Blast CDS nucleotide sequences called by geneMarkS against custom ISKNV and FV3 databases
rule blastn_ISKNV_FV3:
    input:
        rules.geneMarkS.output.fnn
    output:
        blast_ISKNV = join(OUTDIR, "geneMarkS_{samp}/{samp}_GeneMarkS_blastn_ISKNV_AF371960.1.tsv"),
        blast_FV3 = join(OUTDIR, "geneMarkS_{samp}/{samp}_GeneMarkS_blastn_FV3_AY548484.1.tsv"),
    params:
        db_ISKNV = "/labs/kingsley/ambenj/TSIV/analysis/other_virus_sequences/ISKNV_AF371960.1/ISKNV_AF371960.1_orf_nucleotides.fasta",
        db_FV3 = "/labs/kingsley/ambenj/TSIV/analysis/other_virus_sequences/FV3_AY548484.1/FV3_AY548484.1_orf_nucleotides.fasta"
    threads: 4
    resources:
        mem = 16,
        time = 4
    shell: """
        module load ncbi-blast/2.15.0+
        blastn -task dc-megablast \
        -outfmt "6 qseqid sseqid evalue stitle pident length mismatch gapopen qstart qend sstart send bitscore" \
        -query {input} -db {params.db_ISKNV} -num_threads {threads} -out {output.blast_ISKNV} -evalue 0.05

        blastn -task dc-megablast \
        -outfmt "6 qseqid sseqid evalue stitle pident length mismatch gapopen qstart qend sstart send bitscore" \
        -query {input} -db {params.db_FV3} -num_threads {threads} -out {output.blast_FV3} -evalue 0.05
    """

# Blast translated CDS called by geneMarkS against custom ISKNV and FV3 databases
rule blastp_ISKNV_FV3:
    input:
        rules.geneMarkS.output.faa
    output:
        blast_ISKNV = join(OUTDIR, "geneMarkS_{samp}/{samp}_GeneMarkS_blastp_ISKNV_AF371960.1.tsv"),
        blast_FV3 = join(OUTDIR, "geneMarkS_{samp}/{samp}_GeneMarkS_blastp_FV3_AY548484.1.tsv"),
    params:
        db_ISKNV = "/labs/kingsley/ambenj/TSIV/analysis/other_virus_sequences/ISKNV_AF371960.1/ISKNV_AF371960.1_orf_proteins.fasta",
        db_FV3 = "/labs/kingsley/ambenj/TSIV/analysis/other_virus_sequences/FV3_AY548484.1/FV3_AY548484.1_orf_proteins.fasta"
    threads: 4
    resources:
        mem = 16,
        time = 4
    shell: """
        module load ncbi-blast/2.15.0+
        blastp -query {input} -db {params.db_ISKNV} -num_threads {threads} -out {output.blast_ISKNV} -evalue 0.05 \
        -outfmt "6 qseqid sseqid evalue stitle pident length mismatch gapopen qstart qend sstart send bitscore"

        blastp -query {input} -db {params.db_FV3} -num_threads {threads} -out {output.blast_FV3} -evalue 0.05 \
        -outfmt "6 qseqid sseqid evalue stitle pident length mismatch gapopen qstart qend sstart send bitscore"

    """

