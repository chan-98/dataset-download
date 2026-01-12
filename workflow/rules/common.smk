from os.path import join, dirname, basename
from snakemake.utils import validate
from snakemake.exceptions import WorkflowError
from snakemake.io import load_configfile
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import yaml
from typing import List, Dict
from re import match, IGNORECASE
from sys import stderr
import json
from pandas import read_csv


def validate_config(conf: Dict[str, Dict], schemas_path: str) -> None:
    ''' This function validates the config.

    If the config doesn't conform to the schema this function will raise an error.

    Parameters:
        conf: config dictionary to validate
        schemas_path: path to the schemas directory

    Returns:
        None
    '''
    validate(
        conf,
        schema=join(schemas_path, "general.yml")
    )
    return(None)

def get_default_config(conf_path: str) -> Dict:
    return load_configfile(join(conf_path, "default_config.yml"))

def load_filters(source: str, conf_path: str) -> Dict:
    return load_configfile(join(conf_path, f"{source}_filters.yml"))

def get_ensembl_releases(source_subtype: str, config_path: str) -> Dict[str, Dict[str, str]]:
    '''Get Ensembl release information.
    Returns:
        A dictionary containing the different URL for Ensembl archives.
    '''
    with open(os.path.join(config_path, source_subtype+"_releases.yml")) as f:
        # may throw a yaml.YAMLError
        releases = yaml.safe_load(f)
    return releases

def get_allowed_references(source: str, config_path: str) -> Dict[str, Dict[str, str]]:
    '''Get the allowed reference data sets that can be downloaded.
    This comprises:
     - species that can be downloaded from ensembl
     - spike-ins data from vendors

    Parameters:
        source: The source to get the available reference data from
        config_path: Path to the config directory of the pipeline

    Returns:
        A dictionary containing the different datasets available for download for that
        source and their attributes.
    '''
    if source == "ensembl":
        with open(os.path.join(config_path, "ensembl_species.yml")) as f:
            # May raise a yaml.YAMLError
            dat = yaml.safe_load(f)
            return dat
    elif source == "spike-ins":
        with open(os.path.join(config_path, "spike-ins.yml")) as f:
            # May raise a yaml.YAMLError
            dat = yaml.safe_load(f)
            return dat
    else:
        raise WorkflowError(f"Source '{source}' is not available")

def get_selected_species(selected_species: str,
                      allowed_references: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    matched = ""
    for s in allowed_references:
        if match(s, selected_species, IGNORECASE):
            matched = s
            break
    if len(matched) > 0:
        return s
    else:
        raise WorkflowError(f"'{selected_species}' is not in the list of allowed species for {config['source']}")

def get_ensembl_fasta_input(wildcards):
    primary_available = False
    with checkpoints.ensembl_list_files.get(
        source=wildcards.source,
        release=wildcards.release,
        assembly=wildcards.assembly
    ).output[0].open() as f:
        for line in f.readlines():
            l = line.rstrip()
            if l.endswith(f"{wildcards.latin}.{wildcards.assembly}.dna.primary_assembly.fa.gz"):
                primary_available = True
    level = "top_level"
    if primary_available:
        return rules.ensembl_fasta_primary.output
    else:
        return rules.ensembl_fasta_toplevel.output

def get_ensembl_fasta(selected_species: str,
                      allowed_references: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    '''Get ensembl fasta file

    Parameters:
        selected_species: The name of the species to download the reference for
        allowed_references: Dictionary of allowed reference species for ensembl

    Returns:
        The path of the desired fata file
    '''
    assembly = f"{allowed_references[selected_species]['assembly_prefix']}{source_conf['assembly_version']}"
    fasta_path = join(
        outdir, source_conf["species"], config["source"],
        f"release-{source_conf['release_version']}", assembly, "unfiltered",
        f"{allowed_references[selected_species]['latin']}.{assembly}.fa.README")
    return fasta_path

def get_ensembl_gtf(selected_species: str,
                      allowed_references: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    flt = "unfiltered"
    filters = ""
    if source_conf["filtered"]:
        filters = "."
        flt = "filtered"
        for f in source_conf["filters"]:
            key = list(f.keys())[0]
            filters = f"{filters}flt_{key}:{','.join(f[key])}"
    assembly = f"{allowed_references[selected_species]['assembly_prefix']}{source_conf['assembly_version']}"
    gtf_path = join(
        outdir, source_conf["species"], config["source"],
        f"release-{source_conf['release_version']}", assembly, flt,
        f"{allowed_references[selected_species]['latin']}.{assembly}.{source_conf['release_version']}{filters}.gtf.gz")
    return gtf_path


def get_one2one_orthologs_output() -> str:
    return expand(rules.get_one2one_orthologs.output,
        source=source,
        source_species=source_species,
        source_latin=source_latin,
        source_assembly=source_assembly,
        target_species=target_species,
        target_latin=target_latin,
        target_assembly=target_assembly,
        release=release
    )

def get_gene_annotation_sources(config_path: str) -> Dict[str, Dict[str, str]]:
    """
    This function will open the gene annotation source configuration file to 
    extract all the data sources
    """
    with open(os.path.join(config_path, "gene_annotation_sources.yml")) as f:
        # May throw a yaml.YAMLError
        dat = yaml.safe_load(f)
    return dat


def download_gene_annotation_input(wildcards: snakemake.io.Wildcards) -> str:
    category = wildcards.category
    source = wildcards.source
    d = gene_annotation_sources[category]
    ds = d[source]
    return(ds['url'])


def download_gene_annotation_output(wildcards: snakemake.io.Wildcards)  -> str:
    category = wildcards.category
    source = wildcards.source
    d = gene_annotation_sources[category]
    ds = d[source]
    return(ds['filename'])

def download_gene_annotation_description(wildcards: snakemake.io.Wildcards) -> str:
    category = wildcards.category
    source = wildcards.source
    d = gene_annotation_sources[category]
    ds = d[source]
    return(ds['description'])

def get_all_gene_annotations_output(wildcards: snakemake.io.Wildcards) -> List[str]:
    output_files = list()
    for category, v in gene_annotation_sources.items():
        for source, ds in v.items():
            output_files.append(os.path.join(outdir, "genesets", category, source, ds['filename']))
    return output_files

def get_ensembl_biomart(selected_species: str,
                      allowed_references: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    '''Get ensembl BioMart annotations

    Parameters:
        selected_species: The name of the species to download the reference for
        allowed_references: Dictionary of allowed reference species for ensembl

    Returns:
        The path of the desired fata file
    '''
    assembly = f"{allowed_references[selected_species]['assembly_prefix']}{source_conf['assembly_version']}"
    biomart_path = join(
        outdir, source_conf["species"], config["source"],
        f"release-{source_conf['release_version']}", assembly, "unfiltered",
        f"{allowed_references[selected_species]['latin']}.{assembly}.{source_conf['release_version']}.biomart.txt")
    return biomart_path

def get_spike_ins_url(requested_kit: str, allowed_references: Dict[str,str]) -> str:
    if requested_kit in allowed_references.keys():
        kit_data = allowed_references[requested_kit]
        return kit_data['url']
    else:
        kit_list = ', '.join(allowed_references.keys())
        raise WorkflowError(f"'{requested_kit}' is not in the list of allowed kits ({kit_list}) for {config['source']}")

def get_download_spike_ins_output() -> str:
    return expand(rules.download_spike_ins.output, kit=kit)

def get_fastq_screen_build_output(selected_species: str) -> List[str]:
    assembly = f"{allowed_references[selected_species]['assembly_prefix']}{source_conf['assembly_version']}"
    latin = allowed_references[selected_species]['latin']
    bowtie_index_files = expand(rules.fastq_screen_build.output.bt2l_files, 
        species=selected_species,
        source=config["source"],
        release=source_conf['release_version'],
        latin=latin, 
        assembly=assembly)
    log_file = expand(rules.fastq_screen_build.output.log_file,
        species=selected_species,
        source=config["source"],
        release=source_conf['release_version'],
        latin=latin, 
        assembly=assembly)
    bowtie_index_files.extend(log_file)
    return (bowtie_index_files)

def read_sra_source_file(to_read: str) -> List[str]:
    """Reads a file containing a list of SRA accessions to download.

    The file should be formatted as a list, one accessions per row.

    Parameters:
        to_read: File to get the SRA accession list from

    Returns:
        A list of SRA accession.
    """
    # Read accession list
    sra_file = read_csv(to_read, sep="\t", header=None)
    # Return list of accessions after removing white spaces before and after
    return [x.strip() for x in sra_file[0].values.tolist()]

def format_runinfo_input(conf: Dict) -> List[str]:
    """Returns the input files needed by format_runinfo rule.

    Parameters:
        None

    Returns:
        A list of files.
    """
    if len(conf["accessions_to_skip"]) > 0:
        return rules.filter_runinfo.output.filt_meta
    else:
        return rules.get_accession_runinfo.output.unfilt_meta

def format_runinfo_input_runtable_meta(conf: Dict) -> List[str]:
    """Returns the input files needed by format_runinfo rule.

    Parameters:
        None

    Returns:
        A list of files.
    """
    if len(conf["accessions_to_skip"]) > 0:
        return rules.filter_runinfo.output.filt_meta1
    else:
        return rules.get_accession_runinfo.output.unfilt_meta1

def get_pigz_fastq_input(wildcards) -> List[str]:
    """Returns the inputs needed by pigz_fastq.

    The function reads the output from checkpoint vdb_dump to know if the data is FASTQ
    or BAM.

    Parameters:
        None
    
    Returns:
        A list of files.
    """
    with checkpoints.vdb_dump.get(SRR=wildcards.SRR).output[0].open() as f:
        fmt = ""
        for line in f.readlines():
            fmt_match = re.match(r"^FMT\s+\:\s+([a-zA-Z]+)", line)
            if fmt_match is not None:
                fmt = fmt_match.group(1).upper()
        rule_output_dir = ""
        if fmt == "BAM":
            rule_output_dir = "bam2fq"
        elif fmt == "FASTQ":
            rule_output_dir = "fastq_samples"
        else:
            raise WorkflowError(f"Unrecognized format {fmt}")
        return(expand(join(outdir, rule_output_dir, "{SRR}", "{SRR}{ext}"),
                      SRR=wildcards.SRR,
                      ext=wildcards.ext))

def get_bulk_fastq_data(wildcards) -> List[str]:
    """Returns a list of fastq files required for bulk data.

    Reads the runinfo file to know if the output is paired-ended or single-ended

    Parameters:
        None

    Returns:
        A list of files
    """
    # Read runinfo file
    runinfo_file = checkpoints.format_runinfo.get().output.formatted_meta
    runinfo_df = read_csv(runinfo_file, sep="\t")

    # List of files to return
    input_files = list()
    file_extensions = dict()

    # Loop over all the fastq files generated
    for index, row in runinfo_df.iterrows():
        # Remove the path and keep only the filename
        fq = basename(row[0])

        # Get the strandness
        sra_conf[fq] = row["LibraryLayout"]
        file_extensions[fq] = ".fastq"
        if row["LibraryLayout"] == "PAIRED":
            file_extensions[fq] = ["_1.fastq", "_2.fastq"]
    
    # Loop over all the SRR accessions and expand the move_fastq_files rule output accordingly
    for srr_accession in file_extensions:
        input_files.extend(expand(rules.move_fastq_files.output[0], SRR=srr_accession, ext=file_extensions[srr_accession]))
    return input_files

def get_sc_fastq_data(wildcards) -> List[str]:
    """Returns a list of fastq files required for single-cell data.

    Same as get_bulk_fastq_data except the output is expected to always be paired-ended.

    Parameters:
        None

    Returns:
        A list of files
    """
    # Read runinfo file
    runinfo_file = checkpoints.format_runinfo.get().output.formatted_meta
    runinfo_df = read_csv(runinfo_file, sep="\t")

    # List of files to return
    input_files = list()
    file_extensions = dict()

    # Loop over all the fastq files generated
    for index, row in runinfo_df.iterrows():
        # Remove the path and keep only the filename
        fq = basename(row[0])

        # Get the strandness
        sra_conf[fq] = row["LibraryLayout"]
        file_extensions[fq] = ["_1.fastq", "_2.fastq"]
    
    # Loop over all the SRR accessions and expand the move_fastq_files rule output accordingly
    for srr_accession in file_extensions:
        input_files.extend(expand(rules.pigz_fastq.output[0], SRR=srr_accession, ext=file_extensions[srr_accession]))
    return input_files
