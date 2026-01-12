from os.path import join
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule download_spike_ins:
    input:
        HTTP.remote(spike_ins_url, keep_local=True)
    output:
        directory(join(outdir, "spike-ins", kit))
    shell:
        """
        mkdir -p {output} && unzip {input} -d {output}
        """

rule get_spike_ins_fasta_and_gtf:
    input:
        get_download_spike_ins_output()
