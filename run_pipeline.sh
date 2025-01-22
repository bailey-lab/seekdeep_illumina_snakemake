snakemake -s setup_run.smk --cores 8
snakemake -s run_extractor.smk --cores 8
snakemake -s finish_process.smk --cores 8
