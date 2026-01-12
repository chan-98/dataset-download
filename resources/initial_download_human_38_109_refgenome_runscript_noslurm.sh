#!/bin/bash

PIPE='/softwares/grepbio/multiBoard/data-download-pipeline/'
CONFIGFILE='/data/data_download_pipeline/internal/initial_download_human_38_109_refgenome/user_config.yml'

snakemake get_fasta_and_gtf_and_biomart --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile --cores 8 --use-conda
snakemake build_fastq_screen_reference --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile --cores 8 --use-conda