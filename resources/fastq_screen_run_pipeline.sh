#!/bin/bash

PIPE='/softwares/grepbio/multiBoard/data-download-pipeline/'

while read spec; do

        CONFIGFILE='/data/data_download_pipeline/internal/fastq_screen_genomes/'${spec}'_user_config.yml'
        echo "$spec"
        snakemake get_fasta --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile --cores 8 --use-conda
        snakemake get_gtf --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile --cores 8 --use-conda
        snakemake build_fastq_screen_reference --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile --cores 8 --use-conda

done < fastq_screen_genomes_list.txt