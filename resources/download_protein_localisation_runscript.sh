#!/bin/bash

PIPE='/softwares/grepbio/multiBoard/data-download-pipeline'
CONFIGFILE='/data/data_download_pipeline/internal/download_protein_localisation/user_config.yml'

snakemake download_all_annotations --configfile ${CONFIGFILE} --snakefile ${PIPE}/workflow/Snakefile --cores 1