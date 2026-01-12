#!/bin/bash
SBATCH --job-name=data_download      # job name (shows up in the queue)
SBATCH --output=download.log
SBATCH --time=07:10:00         # Walltime (HH:MM:SS)
SBATCH --ntasks=1             # number of tasks (e.g. MPI)
SBATCH --cpus-per-task=1       # number of cores per task (e.g. OpenMP)

PIPE='/home/chandini.v/data-download-pipeline'
CONFIGFILE='/home/chandini.v/data-download-pipeline/user_config.yml'

̈srun snakemake get_fasta --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# snakemake get_fasta --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile ## doesn't work, needs SLURM
# srun snakemake get_gtf --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# srun snakemake get_biomart --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# srun snakemake get_fasta_and_gtf_and_biomart --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# srun snakemake get_spike_ins_fasta_and_gtf --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# srun snakemake download_all_annotations--profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile
# srun snakemake build_fastq_screen_reference --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile $@
# srun snakemake get_fastq_data_from_sra --profile ${PIPE}/slurm --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile $@