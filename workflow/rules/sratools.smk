from os.path import join, dirname, basename

wildcard_constraints:
    SRR="([A-Z0-9]+)",
    ext=r"((_[12])?\.fastq)"

rule get_accession_runinfo:
    output:
        unfilt_meta = temp(join(outdir, "runinfo", "unfiltered_meta.csv")),
        unfilt_meta1 = temp(join(outdir, "runinfo", "unfiltered_meta1.csv"))
    params:
        db = sra_conf["database"],
        query = " OR ".join(accession_list)
    conda:
        join(envs_path, "sratools.yml")
    shell:
        """
        esearch -db {params.db} -query {params.query:q} | \
        efetch -format runinfo > {output.unfilt_meta} &&
        pysradb metadata {params.query:q} --detailed > {output.unfilt_meta1}
        """
 
rule filter_runinfo:
    input:
        unfilt_meta = rules.get_accession_runinfo.output.unfilt_meta,
        unfilt_meta1 = rules.get_accession_runinfo.output.unfilt_meta1
    output:
        filt_meta = temp(join(outdir, "runinfo", "filtered_meta.csv")),
        filt_meta1 = temp(join(outdir, "runinfo", "filtered_meta1.csv"))
    params:
        to_remove = f"^({'|'.join(sra_conf['accessions_to_skip'])})"
    conda:
        join(envs_path, "sratools.yml")
    shell:
        """
        egrep -v {params.to_remove:q} {input.unfilt_meta} > {output.filt_meta} &&
        egrep -v {params.to_remove:q} {input.unfilt_meta1} > {output.filt_meta1}
        """

checkpoint format_runinfo:
    input:
        unformatted_meta = format_runinfo_input(sra_conf),
        unformatted_meta1 = format_runinfo_input_runtable_meta(sra_conf)
    output:
        formatted_meta = join(outdir, "runinfo", "metadata.tsv"),
        formatted_meta1 = join(outdir, "runinfo", "metadata1.tsv")
    shell:
        """
        cat {input.unformatted_meta} | tr ',' '\t'  > {output.formatted_meta} && \
        cat {input.unformatted_meta1} | tr ',' '\t'  > {output.formatted_meta1}
        """

checkpoint vdb_dump:
    output: join(outdir, "SRA_files", "{SRR}", "{SRR}.info")
    conda:
        join(envs_path, "sratools.yml")
    shell:
        """
        vdb-dump --info {wildcards.SRR} --output-file {output}
        """

rule prefetch_sra:
    output:
        sra = temp(join(outdir, "SRA_files", "{SRR}", "{SRR}.sra"))
    params:
        srr_dir = lambda wildcards, output: f"{dirname(dirname(output.sra))}",
        srr_id = lambda wildcards: f"{wildcards.SRR}",
        max_size = sra_conf["prefetch"]["max_size"]
    conda:
        join(envs_path, "sratools.yml")
    shell:
        """
        prefetch {params.srr_id} \
          --max-size {params.max_size} \
          --output-directory {params.srr_dir}
        """

rule fasterq_dump_se:
    input: rules.prefetch_sra.output.sra
    output:
        fq_dir = temp(join(outdir, "fastq_samples", "{SRR}", "{SRR}.fastq"))
    params:
        work_dir = lambda wildcards, input: f"{dirname(dirname(input[0]))}",
        outputdir = lambda wildcards, output: f"{dirname(output.fq_dir)}",
        srr_id = lambda wildcards: f"{wildcards.SRR}",
        temp_dir = sra_conf["fasterq_dump"]["temp"]
    threads: sra_conf["fasterq_dump"]["threads"]
    conda:
        join(envs_path, "sratools.yml")
    shell:
        """
        cd {params.work_dir} && \
        fasterq-dump {params.srr_id} \
            --threads {threads} \
            --temp {params.temp_dir} \
            --outdir {params.outputdir}
        """

rule fasterq_dump_pe:
    input: rules.prefetch_sra.output.sra
    output:
        fq_1 = temp(join(outdir, "fastq_samples", "{SRR}", "{SRR}_1.fastq")),
        fq_2 = temp(join(outdir, "fastq_samples", "{SRR}", "{SRR}_2.fastq"))
    params:
        work_dir = lambda wildcards, input: f"{dirname(dirname(input[0]))}",
        outputdir = lambda wildcards, output: f"{dirname(output.fq_1)}",
        srr_id = lambda wildcards: f"{wildcards.SRR}",
        temp_dir = sra_conf["fasterq_dump"]["temp"]
    threads: sra_conf["fasterq_dump"]["threads"]
    conda:
        join(envs_path, "sratools.yml")
    shell:
        """
        cd {params.work_dir} && \
        fasterq-dump {params.srr_id} \
            --threads {threads} \
            --temp {params.temp_dir} \
            --split-3 \
            --outdir {params.outputdir}
        """

rule sam_dump:
    input: rules.prefetch_sra.output.sra
    output: temp(join(outdir, "fastq_samples", "{SRR}", "{SRR}.bam"))
    threads: 2
    conda:
        join(envs_path, "sratools.yml")
    shell:
        """
        sam-dump --unaligned --header {input} | \
        samtools view -bu - > {output}
        """

rule bam2fq_se:
    input: rules.sam_dump.output[0]
    output: temp(join(outdir, "bam2fq", "{SRR}", "{SRR}.fastq"))
    conda:
        join(envs_path, "samtools.yml")
    threads: sra_conf["bam2fq"]["threads"] - 1
    shell:
        """
        samtools bam2fq \
            --threads {threads} \
            -1 {output} \
            {input}
        """

rule bam2fq_pe:
    input: rules.sam_dump.output[0]
    output:
        fq_1 = temp(join(outdir, "bam2fq", "{SRR}", "{SRR}_1.fastq")),
        fq_2 = temp(join(outdir, "bam2fq", "{SRR}", "{SRR}_2.fastq"))
    conda:
        join(envs_path, "sratools.yml")
    threads: sra_conf["bam2fq"]["threads"] - 1
    shell:
        """
        samtools bam2fq \
            --threads {threads} \
            -1 {output.fq_1} \
            -2 {output.fq_2} \
            {input}
        """

rule pigz_fastq:
    input:
        info = join(outdir, "SRA_files", "{SRR}", "{SRR}.info"),
        to_compress = get_pigz_fastq_input
    output: join(outdir, "fastq_samples", "{SRR}", "{SRR}{ext}.gz")
    threads: sra_conf["pigz"]["threads"]
    conda:
        join(envs_path, "sratools.yml")
    shell:
        """
        pigz -p {threads} -c {input.to_compress} > {output}
        """

rule move_fastq_files:
    input: rules.pigz_fastq.output[0]
    output: join(outdir, "{SRR}{ext}.gz")
    shell:
        """
        mv {input} {output}
        """

# This rule will follow this snakemake logic:
# - If the data is bulk it will be moved to a 'fastq' directory
# - Compress each fastq file separatelly
# - Check if SRA data is BAM of FASTQ format to decide how to procede
#   - If data is BAM
#     - Extract FASTQ files from BAM
#     - Extract BAM from SRA
#   - Else
#     - Extract FASTQ files from SRA
# - Download SRA file using prefetch
# - Filter runinfo file removing accessions_to_skip
# - Download runinfo for all accessioins
if sra_conf["data_type"] == "bulk":
    rule get_fastq_data_from_sra:
        input: get_bulk_fastq_data
else:
    rule get_fastq_data_from_sra:
        input: get_sc_fastq_data
