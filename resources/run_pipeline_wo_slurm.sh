#!/bin/bash

PIPE='/home/minkyung.sung/git/data-download-pipeline'
CONFIGFILE='/home/minkyung.sung/git/data-download-pipeline/resources/user_config.yml'

snakemake get_fasta --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile --use-conda --cores 1
# srun snakemake get_gtf --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# srun snakemake get_biomart --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# srun snakemake get_fasta_and_gtf_and_biomart --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# srun snakemake get_spike_ins_fasta_and_gtf --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# srun snakemake download_all_annotations--profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# srun snakemake build_fastq_screen_reference --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile $@
# srun snakemake get_fastq_data_from_sra --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile $@
