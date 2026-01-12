from os.path import join, splitext
rule fastq_screen_build:
    input:
        join(outdir, "{species}", "{source}", "release-{release}",
            "{assembly}", "unfiltered", "{latin}.{assembly}.dna.toplevel.fa.gz")
    output:
        bt2l_files = multiext(join(outdir,  "indexed_genomes",  "{species}", "{source}", "release-{release}", "{assembly}", "bowtie2", "{latin}.{assembly}"), ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"),
        log_file = join(outdir, "indexed_genomes", "{species}", "{source}", "release-{release}", "{assembly}", "bowtie2", "{latin}.{assembly}_README.txt")

    message: "Using: {wildcards.latin}.{wildcards.assembly}.{wildcards.release}.dna.toplevel.fa.gz"
    params:
        bt2l_index_prefix = lambda wildcards, output: splitext(splitext(output.bt2l_files[0])[0])[0],
        bt2l_threads = config['bowtie2_params']['num_threads'] - 2 if config['bowtie2_params']['num_threads'] > 2 else config['bowtie2_params']['num_threads']
    resources: 
        mem=config['bowtie2_params']['mem']
    threads:
        config['bowtie2_params']['num_threads']
    conda:
        join(envs_path, "fastq_screen_build.yml")
    shell:
        """
        bowtie2-build --large-index --threads={threads} '{input}' '{params.bt2l_index_prefix}' && \
        echo "Generated index on {wildcards.latin}.{wildcards.assembly} from {wildcards.source} using fastq_screen_build rule of Download data pipeline on "$(date) > {output.log_file}
        """

rule build_fastq_screen_reference:
    input:
        get_fastq_screen_build_output(selected)
