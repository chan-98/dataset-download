from os.path import join, dirname

rule get_one2one_orthologs:
    input:
        source_biomart = join(outdir, "{source_species}", "{source}", "release-{release}", "{source_assembly}", "unfiltered", "{source_latin}.{source_assembly}.{release}.biomart.txt"),
        target_biomart = join(outdir, "{target_species}", "{source}", "release-{release}", "{target_assembly}", "unfiltered", "{target_latin}.{target_assembly}.{release}.biomart.txt")
    output:
        join(outdir, "{source_species}", "{source}", "release-{release}", "{source_assembly}", "orthologs", "{target_species}", "{target_latin}.{target_assembly}.{source_latin}.{source_assembly}_one2one_genemap.txt")
    message:
        "Get {source_species} {target_species} orthologs from Biomart for {source} {release}"
    params:
        source_name = lambda wildcards: f"{source_species}".lower(),
        target_name = lambda wildcards: f"{target_species}".lower(),
        release = lambda wildcards: f'{release}',
        biomart_url = lambda wildcards: f'{biomart_url}',
        basedir = workflow.basedir
    conda:
        join(envs_path, "etl.yml")
    shell:
        """
        Rscript {params.basedir}/scripts/etl/get_one2one_orthologs.R \
                --source_genes {input.source_biomart} \
                --target_genes {input.target_biomart} \
                --source_latin_prefix {source_latin_prefix} \
                --target_latin_prefix {target_latin_prefix} \
                --source_name {params.source_name} \
                --target_name {params.target_name} \
                --biomart {params.biomart_url} \
                --release_version {params.release} \
                --output {output}
        """

rule get_orthologs:
    input:
        get_one2one_orthologs_output()