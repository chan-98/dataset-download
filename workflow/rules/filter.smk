from os.path import join
import re


rule filter_gtf:
    input:
        join(outdir, "{species}", "{source}", "release-{release}", "{assembly}",
            "unfiltered", "{latin}.{assembly}.{release}.gtf")
    output:
        join(outdir, "{species}", "{source}", "release-{release}", "{assembly}",
            "filtered", "{latin}.{assembly}.{release}.{filters}.gtf")
    params:
        script_path = join(scripts_path, "filter_gtf.py"),
        filters = "{filters}"
    shell:
        """
        python {params.script_path} \
            --input '{input}' \
            --output '{output}' \
            --filters '{params.filters}'
        """

rule gunzip_gtf:
    input: "{file_path}.gtf.gz"
    output: "{file_path}.gtf"
    shell:
        """
        gzip -dc '{input}' > '{output}'
        """

rule gzip_gtf:
    input: rules.filter_gtf.output
    output: f"{rules.filter_gtf.output[0]}.gz"
    shell:
        """
        gzip -c '{input}' > '{output}'
        """