# Data download pipeline

This Snakemake pipeline will download data and store it systematicaly.

**Features:**

- Download genomes (FASTA and GTF) from ensembl
- Filter GTf file
- Extract gene annotations for downstream operations

**Available rules:**

- get_fasta
- get_gtf
- get_biomart
- get_fasta_and_gtf_and_biomart
- build_fastq_screen_reference
- get_fastq_data_from_sra

**Available sources:**

- ensembl
- sra

**Available species**

- ensembl
  - Human
  - Mouse
  - Macaque
  - Crab-eating macaque

## Usage

- [Reference genomes](docs/reference_genomes.md)
- [Gene and cell annotations](docs/gene_annotations.md)

To use this pipeline you just need to copy the files in the `resources/` directory.

### user_config.yml

In the `user_config.yml` file you can specify the source of the data
to download

#### Reference genome

If the source is set to `ensembl`, you can then specify what assembly version and
emsembl release version to download.

Word of caution: the pipeline contains a configuration file called
`ensembl_releases.yml` that points each release of Ensembl to the corresponding
Ensembl archive URL except for the latest that will point to 'www.ensembl.org'.
It is therefore imperative to update this file on a regular basis to avoid pointing
a release number to the wrong site. You can suscribe to the Ensembl mailing list
to get the latest news (<announce-join@ensembl.org>).

If `filtered` is set to `True` the resulting filtered GTF will only contain lines
matching the filters listed.

To specify filters you can add a `filters` section as shown in the example below:

```yml
ensembl_params:
  filtered: False
  # Filters to apply to the GTF file (lines matching these will be kept)
  filters:
    - gene_biotype:
      - "+protein_coding"
    - transcript_biotype:
      - "-lncRNA"
```

In this example the pipeline will output a file containing the only `protein_coding`
genes and non `lncRNA` transcripts.

#### Barnyard genome

The Barnyard genome is a concatenation of Human and Mouse genome annotations.
Human and mouse reference dataset, GRCh38 and mm10 can be downloaded directly from the
[10x genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38mm10_2020A).
Chromosome assembly sequences and gene annotations in gtf format will be provided by default.
To generate the BioMart annotation file that is required by our pipelines, you can run the
following command (currently not wrapped as a pipeline rule):

```{shell}
python ./data-download-pipeline/workflow/scripts/create_barnyard_biomart.py --human=Homo_sapiens.GRCh38.98.biomart.txt --mouse=Mus_musculus.GRCm38.98.biomart.txt -o Barnyard.biomart.txt
```

#### Gene and cell annotations

To download reference gene annotation data, you need to create a configuration file
as simple as:

```
source: "other"
destination: "/data/public/4.Reference_genomes/references"
create_gene_sets: True
```

The source is set to 'other' because there are many sources for gene and cell annotations.
Details about each source is available in [Gene and cell annotations](docs/gene_annotations.md).

'create_gene_sets' is set to true to indicate to the pipeline that it needs to
include the rules that will update the gene and cell annotations.

The following rules can be executed with Snakemake:

- download_all_annotations: for any source defined in [config/gene_annotation_sources.yml](config/gene_annotation_sources.yml)
- get_homo_sapiens_surface_gene_list: to extract all the cell surface markers from one of the sources above
- get_mt_rp_gene_list: to get Human MT gene and Ribosomal protein list using Biomart and GTF gene annotations
- get_homo_sapiens_uniprot_subcellular_location: to retrieve all subcellular location information from UniProt

### One to One ortholog mapping

The single-cell pipeline rely on the availability of Mouse to Human
1 to 1 gene ortholog mappings. To generate mappings between 2 species for a
defined release of Ensembl, you need to prepare a configuration file with the
following properties:

```
source: "other"
destination: "/data/public/4.Reference_genomes/references"
orthologs:
    source_species: Human
    source_assembly_version: 38
    target_species: Mouse
    target_assembly_version: 38
    release_version: 98
```

The corresponding rule to execute with Snakemake is:

```
srun snakemake get_orthologs ...
```

In the above example, results will be stored in:

```
/data/public/4.Reference_genomes/references/Human/ensembl/release-98/GRCh38/orthologs/Mouse/Mus_musculus.GRCm38.Homo_sapiens.GRCh38_one2one_genemap.txt
```

Please note that in its current implementation, you need to run the pipeline to
download the required source and target species genome information as specified
in the Reference genome section above.

#### Spike-ins sequences

Variation in RNA expression data can be attributed to a variety of factors including the quality of the starting material, the level of cellularity and RNA yield, the platform employed, and the person performing the experiment. To control for these sources of variability, a common set of external RNA controls has been developed by the External RNA Controls Consortium (ERCC), an ad-hoc group of academic, private, and public organizations hosted by the National Institute of Standards and Technology (NIST). The controls consist of a set of unlabelled, polyadenylated transcripts designed to be added to an RNA analysis experiment after sample isolation, in order to measure against defined performance criteria.
In RNA-Seq analyses, adding pre-determined quantity of synthetic RNA sequences (spike-ins) to samples is a popular way to verify the experimental pipeline, determine quantification accuracy and for normalisation of differential expression.

The most commonly used spike-ins are the ERCC spike-ins.
Different versions of the sequences are available. For instance, the ERCC Spike-In Control Mixes are commercially available, pre-formulated blends of 92 transcripts, derived and traceable from NIST-certified DNA plasmids. The transcripts are designed to be 250 to 2,000 nt in length, which mimic natural eukaryotic mRNAs.

If the source is set to `spike-ins`, you can then specify what spike-ins kit you
want to download. Currently, we only support the ERCC92 kit from ThermoFisher.

```yml
# Where to download the data from
source: "spike-ins"

# Where to store the data
destination: "/data/public/4.Reference_genomes/references"

# ensembl parameters
spike_ins_params:
  kit: "ERCC92"
```

#### FastQ Screen

If `build_fastq_screen_index` is set to `True`, a bowtie2 index will be created
and placed in the destination directory under `index/bowtie`

If you want to see all possible parameters you can refer to [default_config.yml](./config/default_config.yml).
Copy the parameter you wish to change into your config and change it as desired.

Only species in [ensembl_species.yml](./config/ensembl_species.yml) can be downloaded.

#### SRA FastQ data

If you wish to download data hosted in SRA you will need to set the `source` to 'sra' and the destination to the directory you wish to store the downloaded FastQ files in.

```yml
# Where to download the data from
source: "sra"

# Where to store the data
destination: "/data/ngs/raw/public/053.PRJNA7878"
```

After setting these two parameters, you need to set the parameters in the `sra_params` section.

`data_type` will control if, for the case of 'bulk' data, all the FastQ files will be collated into one single directory called 'fastq'.
Or, if `data_type` is set to 'single-cell' the FastQ files will be left in their sample directory as follows 'fastq_samples/SAMPLE/'.

The `source_type` can be an "accession list" or a "file". In the case of a file you just need to write the path to the file containing the list of accessions in `accessions` and the file should list the accessions as shown below.

```
SRRXXXX1
SRRXXXX2
SRRXXXX3
```

In the case of an an "accession list" you can list the accessions you wish to download in `accessions` as follows:

```yml
# SRA parameters
sra_params:
  # Data type. Can be one of "bulk" or "single-cell"
  data_type: "bulk"
  # Can be one of "file", containing a list of accessions, or "accession list"
  source_type: "accession list"
  # If the source_type is "file" write the path to said file ("/data/path/to/file")
  # If it's "accession list" write the accessions to download (["SRRXX1", "SRRXX2"])
  accessions: ["SRRXX1", "SRRXX2"]
  # This can also be written as a non-pythonic list
  accessions:
    - "SRRXX1"
    - "SRRXX2"
```

Please note that in the accession you can also write project (PRJNA) and experiment accessions (SRX).

In case you want to download all the SRR runs of a project except for a few. You can use the `accessions_to_skip` parameter.

The pipeline will download the runinfo for all the accessions provided in `accessions` and will then filter out the accessions in `accessions_to_skip`.

All fastq files are compressed using `pigz` as a final step.

To run the pipeline you just need to edit the `run_pipeline.sh` file to run the rule `get_fastq_data_from_sra`.

### run-pipeline.sh

In `run_pipeline.sh` you can specify the rule to run (see available rules [Available rules](#data-download-pipeline)).

To run the pipeline you just need to specify which config file to use, activate the
snakemake environment and run the script using sbatch:

```sh
# Activate snakemake environmnet
conda activate snamkemake
# Run script
sbatch run_pipeline.sh
```
