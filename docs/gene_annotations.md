
[[_TOC_]]

# Introduction

The first principle is that any extraction and transformation of existing data
in CI should be traceable and reproducible.
The second principle when it comes to genomic data is that the lineage of the
transformation should indicate (whenever possible):

- The genome build
- The ensembl version or
- The version of the original sources(s) e.g. UniProt release, etc.

For any genome reference, the official source will be the data downloaded by the
data download pipeline.

Information about ETL and data governance is provided in the
[library](http://grepbio.sharepoint.com)

# Data sources

## Surfaceome

The surfaceome can be defined as all plasma membrane proteins that have at least
one amino acid residue exposed to the extracellular space. Therefore, the
surfaceome is a subset of the plasma membrane proteome, which is a subset of the
membrane proteome, the entirety of all membrane proteins
( [ref](https://www.pnas.org/doi/10.1073/pnas.1808790115) ).

All surfaceome sources are defined in `config/gene_annotation_sources.yml`

| Source                        | Description                            | Species |
| ----------------------------- | -------------------------------------- | ------- |
| ms_cell_surface_protein_atlas | CSPA validated surfaceome proteins     | Human   |
| in_silico_human_surfaceome    | in silico surfaceome of 2,886 proteins | Human   |

## Subcellular location

Subcellular location sources are defined in `config/gene_annotation_sources.yml`

| Source                  | Description                                                                  |
| ----------------------- | ---------------------------------------------------------------------------- |
| the_human_protein_atlas | subcellular location of proteins based on immunofluorescently stained cells  |

There is one exception with the UniProt subcellular location because extracting
the information requires to call a REST API. A separate script has been developed
for that purpose as indicated in the ETL section below.

## Genes

The HGNC is the authority for approving human gene symbols and corresponding
descriptive gene names. It is common that other sources of information rely on
human gene symbols. Having a dump of all current gene symbols associated to
Ensembl Ids is therefore desirable.

| Source | Description                                                           |
| ------ | --------------------------------------------------------------------- |
| hgnc   | Human gene symbols and corresponding descriptive gene names from HGNC |

## Cell markers

PanglaoDB contains pre-processed and pre-computed analyses from more than 1054
single-cell experiments covering most major single cell platforms and protocols,
based on more than 4 million cells from a wide range of tissues and organs.
PanglaoDB established a community-curated cell-type marker compendium,
containing more than 6000 gene-cell-type associations, as a resource for
automatic annotation of cell types.

| Source      | Description                                                                                          |
| ----------- | ---------------------------------------------------------------------------------------------------- |
| panglaodb   | community-curated cell-type marker compendium, containing more than 6000 gene-cell-type associations |

## Snakemake rule

| Rule                               | Description                        | Species | Usage      |
| ---------------------------------- | ---------------------------------- | ------- | ---------- |
| download_all_annotations           | Download all annotations           | All     | scPipeline |

# ETL

Transformation rules depend on the extraction of gene annotations from
Ensembl BioMart. These annotations files are downloaded by the pipeline and
made available in:

```
/data/public/4.Reference_genomes/references
```

| Rule                               | Description                        | Species | Usage      |
| ---------------------------------- | ---------------------------------- | ------- | ---------- |
| get_homo_sapiens_surface_gene_list | CSPA validated surfaceome proteins | Human   | scPipeline |
| get_mt_rp_gene_list | Mitochondrial and Ribosomal genes Filtering | Human/Mouse | scPipeline |
| get_homo_sapiens_uniprot_subcellular_location | UniProt subcellular location | Human | scPipeline |
