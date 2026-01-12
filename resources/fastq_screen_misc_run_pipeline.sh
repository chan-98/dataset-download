mamba env create -n fastq_screen -f /softwares/grepbio/multiBoard/data-download-pipeline/workflow/envs/fastq_screen.yml
mamba activate fastq_screen
fastq_screen --get_genomes --outdir /data/public/references/indexed_genomes/other_orgs