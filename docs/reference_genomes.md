**Available rules:**

- get_fasta
- get_gtf
- get_fasta_and_gtf

**Available sources:**

- ensembl

**Available species**

- ensembl
  - Human
  - Mouse
  - Macaque
  - Crab-eating macaque

### user_config.yml

In the `user_config.yml` file you can specify what assembly version and emsembl release
version to download.

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

If `build_fastq_screen_index` is set to `True`, a bowtie2 index will be created
and placed in the destination directory under `index/bowtie`

If you want to see all possible parameters you can refer to [default_config.yml](./config/default_config.yml).
Copy the parameter you wish to change into your config and change it as desired.

Only species in [ensembl_species.yml](./config/ensembl_species.yml) can be downloaded.

## "Error source is not available" what does it mean?

This means that the way to download this species data hasn't been added to the [ensembl_species.yml](./config/ensembl_species.yml). Only species added to this file can be downloaded using ensembl.
