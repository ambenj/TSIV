# TSIV

## Assembly of TSIV genome
1. Quality check raw reads
2. Trim reads
3. Map reads to stickleback reference genome
4. Extract unmapped reads
5. Assembly reads using SPADES
6. Evaluate SPADES assembly with coverage

```bash
# Run SPADES assembly pipeline
module load snakemake
snakemake --configfile configs/240904_spades_assembly_from_raw_STMY_2012_42.yaml --snakefile scripts/assembly_and_annotation/spades_assembly_from_raw.smk --profile scg --jobs 300 --restart-times 0 --rerun-incomplete
```

