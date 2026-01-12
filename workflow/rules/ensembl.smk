from os.path import join, dirname

checkpoint ensembl_list_files:
    output:
        join(outdir, selected, "{source}", "release-{release}",
            "{assembly}", "available_list.txt")
    message: "Fetch list of available files"
    params:
        release = f"release-{source_conf['release_version']}",
        latin_lower = allowed_references[selected]["latin"].lower(),
        source_subtype = {config['source_subtype']}
    shell:
        """
        if [[ {params.source_subtype} = "ensembl" ]]
        then
            wget 'https://ftp.ensembl.org/pub/{params.release}/fasta/{params.latin_lower}/dna/CHECKSUMS' -O '{output}'
        elif [[ {params.source_subtype} = "ensembl_plants" ]]
        then
            wget 'https://ftp.ebi.ac.uk/ensemblgenomes/pub/{params.release}/plants/fasta/{params.latin_lower}/dna/CHECKSUMS' -O '{output}'
        elif [[ {params.source_subtype} = "ensembl_bacteria" ]]
        then
            wget 'https://ftp.ebi.ac.uk/ensemblgenomes/pub/{params.release}/bacteria/fasta/bacteria_0_collection/{params.latin_lower}/dna/CHECKSUMS' -O '{output}'
        fi

        """


rule ensembl_fasta_primary:
    input:
        join(outdir, selected, "{source}", "release-{release}",
            "{assembly}", "available_list.txt")
    output:
        join(outdir, selected, "{source}", "release-{release}",
            "{assembly}", "unfiltered", "{latin}.{assembly}.dna.primary_assembly.fa.gz")

    message: "Dowloading: {wildcards.latin}.{wildcards.assembly}.{wildcards.release}.dna.primary_assembly.fa.gz"
    params:
        release = f"release-{source_conf['release_version']}",
        latin = allowed_references[selected]["latin"],
        latin_lower = allowed_references[selected]["latin"].lower(),
        assembly = f"{allowed_references[selected]['assembly_prefix']}{source_conf['assembly_version']}",
        source_subtype = {config['source_subtype']}
    conda:
        join(envs_path, "ensembl.yml")
    shell:
        """
        if [[ {params.source_subtype} = "ensembl" ]]
        then
            wget 'http://ftp.ensembl.org/pub/{params.release}/fasta/{params.latin_lower}/dna/{params.latin}.{params.assembly}.dna.primary_assembly.fa.gz'\
            -O '{output}'
        elif [[ {params.source_subtype} = "ensembl_plants" ]] 
        then
            wget 'http://ftp.ebi.ac.uk/ensemblgenomes/pub/{params.release}/plants/fasta/{params.latin_lower}/dna/{params.latin}.{params.assembly}.dna.primary_assembly.fa.gz'\
            -O '{output}'
        elif [[ {params.source_subtype} = "ensembl_bacteria" ]]
        then
            wget 'http://ftp.ebi.ac.uk/ensemblgenomes/pub/{params.release}/bacteria/fasta/bacteria_0_collection/{params.latin_lower}/dna/{params.latin}.{params.assembly}.dna.primary_assembly.fa.gz'\
            -O '{output}'
        fi
        """

rule ensembl_fasta_transcriptome:
    input:
    output:
        join(outdir, selected, "{source}", "release-{release}",
            "{assembly}", "unfiltered", "{latin}.{assembly}.cdna.all.fa.gz")

    message: "Dowloading: {wildcards.latin}.{wildcards.assembly}.{wildcards.release}.cdna.all.fa.gz"
    params:
        release = f"release-{source_conf['release_version']}",
        latin = allowed_references[selected]["latin"],
        latin_lower = allowed_references[selected]["latin"].lower(),
        assembly = f"{allowed_references[selected]['assembly_prefix']}{source_conf['assembly_version']}",
        source_subtype = {config['source_subtype']}
    conda:
        join(envs_path, "ensembl.yml")
    shell:
        """
        if [[ {params.source_subtype} = "ensembl" ]]
        then
            wget 'http://ftp.ensembl.org/pub/{params.release}/fasta/{params.latin_lower}/cdna/{params.latin}.{params.assembly}.cdna.all.fa.gz'\
            -O '{output}'
        elif [[ {params.source_subtype} = "ensembl_plants" ]]
        then
            wget 'http://ftp.ebi.ac.uk/ensemblgenomes/pub/{params.release}/plants/fasta/{params.latin_lower}/cdna/{params.latin}.{params.assembly}.cdna.all.fa.gz'\
            -O '{output}'
        elif [[ {params.source_subtype} = "ensembl_bacteria" ]]
        then
            wget 'http://ftp.ebi.ac.uk/ensemblgenomes/pub/{params.release}/bacteria/fasta/bacteria_0_collection/{params.latin_lower}/cdna/{params.latin}.{params.assembly}.cdna.all.fa.gz'\
            -O '{output}'
        fi
        """

rule ensembl_fasta_toplevel:
    input:
        join(outdir, selected, "{source}", "release-{release}",
            "{assembly}", "available_list.txt")
    output:
        join(outdir, selected, "{source}", "release-{release}",
            "{assembly}", "unfiltered", "{latin}.{assembly}.dna.toplevel.fa.gz")
    message: "Dowloading: {wildcards.latin}.{wildcards.assembly}.{wildcards.release}.dna.toplevel.fa.gz"
    params:
        release = f"release-{source_conf['release_version']}",
        latin = allowed_references[selected]["latin"],
        latin_lower = allowed_references[selected]["latin"].lower(),
        assembly = f"{allowed_references[selected]['assembly_prefix']}{source_conf['assembly_version']}",
        source_subtype = {config['source_subtype']}
    conda:
        join(envs_path, "ensembl.yml")
    shell:
        """
        if [[ {params.source_subtype} = "ensembl" ]]
        then
            wget 'http://ftp.ensembl.org/pub/{params.release}/fasta/{params.latin_lower}/dna/{params.latin}.{params.assembly}.dna.toplevel.fa.gz'\
        -O '{output}'
        elif [[ {params.source_subtype} = "ensembl_plants" ]]
        then
            wget 'http://ftp.ebi.ac.uk/ensemblgenomes/pub/{params.release}/plants/fasta/{params.latin_lower}/dna/{params.latin}.{params.assembly}.dna.toplevel.fa.gz'\
            -O '{output}'
        elif [[ {params.source_subtype} = "ensembl_bacteria" ]]
        then
            wget 'http://ftp.ebi.ac.uk/ensemblgenomes/pub/{params.release}/bacteria/fasta/bacteria_0_collection/{params.latin_lower}/dna/{params.latin}.{params.assembly}.dna.toplevel.fa.gz'\
            -O '{output}'
        fi
        """


rule get_ensembl_fasta:
    input:
        rules.ensembl_fasta_transcriptome.output,
        get_ensembl_fasta_input
    output:
        join(outdir, "{species}", "{source}","release-{release}", "{assembly}",
            "unfiltered", "{latin}.{assembly}.fa.README")
    shell:
        """
        echo '# Created: ' > '{output}' &&
        date >> '{output}'
        """
    

rule get_fasta:
    input:
        get_ensembl_fasta(selected, allowed_references)

rule get_ensembl_biomart:
    output:
        join(outdir, "{species}", "{source}", "release-{release}", "{assembly}",
            "unfiltered", "{latin}.{assembly}.{release}.biomart.txt")
    message:
        "Get Biomart annotations"
    params:
        release = f"{source_conf['release_version']}",
        biomart_url = ensembl_releases[f"{source_conf['release_version']}"],
        biomart_gene_dataset = allowed_references[selected]["biomart_gene_dataset"],
        biomart_external_symbol = allowed_references[selected]["biomart_external_symbol"],
        biomart_external_id = allowed_references[selected]["biomart_external_id"],
        latin = allowed_references[selected]["latin"],
        latin_lower = allowed_references[selected]["latin"].lower(),
        assembly = f"{allowed_references[selected]['assembly_prefix']}{source_conf['assembly_version']}"
    conda:
        join(envs_path, "ensembl.yml")
    shell:
        """
        wget -O '{output}' '{params.biomart_url}/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query> \
  <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" > \
    <Dataset name = "{params.biomart_gene_dataset}" interface = "default" > \
      <Attribute name = "ensembl_gene_id" /> \
      <Attribute name = "ensembl_transcript_id" />\
      <Attribute name = "description" />\
      <Attribute name = "external_gene_name" />\
      <Attribute name = "external_gene_source" />\
      <Attribute name = "gene_biotype" />\
      <Attribute name = "transcript_biotype" />\
      <Attribute name = "entrezgene_id" />\
      <Attribute name = "{params.biomart_external_symbol}" />\
      <Attribute name = "{params.biomart_external_id}" />\
      <Attribute name = "ensembl_gene_id_version" />\
      <Attribute name = "ensembl_transcript_id_version" />\
      <Attribute name = "ensembl_peptide_id" />\
    </Dataset>\
  </Query>'
        sed -i -e 1i'ensembl_gene_id\tensembl_transcript_id\tdescription\tgene_name\texternal_gene_source\tgene_biotype\ttranscript_biotype\tentrezgene_id\t{params.biomart_external_symbol}\t{params.biomart_external_id}\tensembl_gene_id_version\tensembl_transcript_id_version\tensembl_peptide_id' '{output}'
        """
rule get_ensembl_gtf:
    output:
        join(outdir, "{species}", "{source}", "release-{release}", "{assembly}",
            "unfiltered", "{latin}.{assembly}.{release}.gtf.gz")
    params:
        release = f"release-{source_conf['release_version']}",
        release_number = f"{source_conf['release_version']}",
        latin = allowed_references[selected]["latin"],
        latin_lower = allowed_references[selected]["latin"].lower(),
        assembly = f"{allowed_references[selected]['assembly_prefix']}{source_conf['assembly_version']}",
        source_subtype = {config['source_subtype']}
    conda:
        join(envs_path, "ensembl.yml")
    shell:
        """
        if [[ {params.source_subtype} = "ensembl" ]]
        then
            wget 'http://ftp.ensembl.org/pub/{params.release}/gtf/{params.latin_lower}/{params.latin}.{params.assembly}.{params.release_number}.gtf.gz'\
            -O '{output}'
        elif [[ {params.source_subtype} = "ensembl_plants" ]]
        then
            wget 'http://ftp.ebi.ac.uk/ensemblgenomes/pub/{params.release}/plants/gtf/{params.latin_lower}/{params.latin}.{params.assembly}.{params.release_number}.gtf.gz'\
            -O '{output}'
        elif [[ {params.source_subtype} = "ensembl_bacteria" ]]
        then
            wget 'http://ftp.ebi.ac.uk/ensemblgenomes/pub/{params.release}/bacteria/gtf/bacteria_0_collection/{params.latin_lower}/{params.latin}.{params.assembly}.{params.release_number}.gtf.gz'\
            -O '{output}'
        fi
        """

rule get_gtf:
    input:
        get_ensembl_gtf(selected, allowed_references)

rule get_fasta_and_gtf_and_biomart:
    input:
        get_ensembl_fasta(selected, allowed_references),
        get_ensembl_gtf(selected, allowed_references),
        get_ensembl_biomart(selected, allowed_references)
