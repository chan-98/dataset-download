import os

#------------------------------------------
# Homo Sapiens - UniProt subcellular location
#------------------------------------------
rule get_homo_sapiens_uniprot_subcellular_location:
    input:
        hgnc = os.path.join(outdir, "genesets", "genes", "hgnc", "hgnc.tsv"),
    output:
        os.path.join(outdir, "genesets", "etl", "Human", "uniprot", "subcellular_location.tsv")
    message: 
        "Get UniProt Human proteins subcellular location"
    params:
        basedir = workflow.basedir
    conda:
        os.path.join(envs_path, "etl_python.yml")
    shell:
        """
        python {params.basedir}/scripts/etl/get_UniProt_subcellular_location.py \
                --genenames {input.hgnc} \
                --output {output}
        """

#------------------------------------------
# Homo Sapiens - Surfaceome list
#------------------------------------------
rule get_homo_sapiens_surface_gene_list:
    input:
        surfaceome = os.path.join(outdir, "genesets", "surfaceome", "in_silico_human_surfaceome", "table_S3_surfaceome.xlsx"),
        biomart = '/data/public/4.Reference_genomes/references/Human/ensembl/release-98/GRCh38/unfiltered/Homo_sapiens.GRCh38.98.biomart.txt',
        fantom = '/projects/ci/results/data/raw/FFcage52_GRCh38_tpm.txt'
    output:
        os.path.join(outdir, "genesets", "etl", "Human", "release-98", "ff_bm_surface_annot.tsv")
    message: 
        "Get surface gene list using Biomart gene annotation and surfaceome dump"
    params:
        basedir = workflow.basedir
    conda:
        os.path.join(envs_path, "etl_openssl3.yml")
    shell:
        """
        Rscript {params.basedir}/scripts/etl/getSurfaceGeneList.R \
                --biomart {input.biomart} \
                --surfaceome {input.surfaceome} \
                --fantom {input.fantom} \
                --output {output}
        """


#------------------------------------------
# Homo Sapiens - MT and RP genes
#------------------------------------------
rule get_homo_sapiens_mt_rp_gene_list:
    input:
        biomart = '/data/public/4.Reference_genomes/references/Human/ensembl/release-98/GRCh38/unfiltered/Homo_sapiens.GRCh38.98.biomart.txt',
        gtf = '/data/public/4.Reference_genomes/references/Human/ensembl/release-98/GRCh38/unfiltered/Homo_sapiens.GRCh38.98.gtf.gz'
    output:
        os.path.join(outdir, "genesets", "etl", "Human", "release-98","GRCh38", "MT_and_RP_genes_human.txt")
    message: 
        "Get Human MT gene and Ribosomal protein list using Biomart and GTF gene annotations"
    params:
        basedir = workflow.basedir
    conda:
        os.path.join(envs_path, "etl_openssl3.yml")
    shell:
        """
        Rscript {params.basedir}/scripts/etl/get_MT_and_RP_genes.R \
                --biomart {input.biomart} \
                --gtf {input.gtf} \
                --output {output}
        """

#------------------------------------------
# Mus Musculus - MT and RP genes
#------------------------------------------
rule get_mus_musculus_mt_rp_gene_list:
    input:
        #biomart = '/projects/ci/results/ext_annotation/GRCm38_geneannot.txt',
        biomart = '/data/public/4.Reference_genomes/references/Mouse/ensembl/release-98/GRCm38/unfiltered/Mus_musculus.GRCm38.98.biomart.txt',
        gtf = '/data/public/4.Reference_genomes/references/Mouse/ensembl/release-98/GRCm38/unfiltered/Mus_musculus.GRCm38.98.gtf.gz'
    output:
        os.path.join(outdir, "genesets", "etl", "Mouse", "release-98", "GRCm38", "MT_and_RP_genes_mouse.txt")
    message: 
        "Get Mouse MT gene and Ribosomal protein list using Biomart and GTF gene annotations"
    params:
        basedir = workflow.basedir
    conda:
        os.path.join(envs_path, "etl_openssl3.yml")
    shell:
        """
        Rscript {params.basedir}/scripts/etl/get_MT_and_RP_genes.R \
                --biomart {input.biomart} \
                --gtf {input.gtf} \
                --output {output}
        """

rule get_mt_rp_gene_list:
    message:
        "Get all MT gene and Ribosomal protein lists (Hs and Mm)"
    input:
        rules.get_mus_musculus_mt_rp_gene_list.output,
        rules.get_homo_sapiens_mt_rp_gene_list.output

rule download_gene_annotation_file:
    output:
        os.path.join(outdir, "genesets", "{category}", "{source}", "{filename}")
    message:
        "Downloading {wildcards.source} data for '{wildcards.category}' category"
    params:
        download_link = lambda wildcards: download_gene_annotation_input(wildcards),
        description = lambda wildcards: download_gene_annotation_description(wildcards),
        output_dir = os.path.join(outdir, "genesets", "{category}", "{source}"),
        user_agent = "Mozilla/5.0 (Windows NT 6.1; Win64; x64; rv:59.0) Gecko/20100101 Firefox/59.0" 
    conda:
        os.path.join(envs_path, "etl_openssl3.yml")
    shell:
        """
        mkdir -p {params.output_dir} && curl -A "{params.user_agent}" -L '{params.download_link}' --output '{output}' && echo '{params.description}. Downloaded on '$(date) > {params.output_dir}/README.txt
        """

rule download_all_annotations:
    message:
        "Downloading all gene/cell marker annotations"
    input:
        get_all_gene_annotations_output
