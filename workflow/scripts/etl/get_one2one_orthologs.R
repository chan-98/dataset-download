library(biomaRt)
library(optparse)

###############################################################################
#' Get one to one orthologs between 2 species from Biomart
#'
#' @name get_one2one_orthologs.R
#' @author Minkyung Sung, Kalaivani Raju, Gautier Koscielny
#' @date 9 May 2023
#' @description This script takes 2 species (one source one target) and extracts
#' the one to one orthologs from BioMart
#'
###############################################################################

# Automatic string to factor conversion introduces non-reproducibility
options(stringsAsFactors = FALSE)

#' Function to parse arguments
#'
#' @examples
#' args <- get_args()
get_args <- function() {
    option_list <- list(
        make_option(
            "--source_genes",
            default = "",
            type = "character",
            help = "Input source species biomart filename"
        ),
        make_option(
            "--target_genes",
            default = "",
            type = "character",
            help = "Input target species biomart filename"
        ),
        make_option(
            "--source_latin_prefix",
            default = "",
            type = "character",
            help = "Input source latin name (e.g. hsapiens)"
        ),
        make_option(
            "--target_latin_prefix",
            default = "",
            type = "character",
            help = "Input target latin name (e.g. mmusculus)"
        ),
        make_option(
            "--source_name",
            default = "",
            type = "character",
            help = "Input source species name (e.g. human)"
        ),
        make_option(
            "--target_name",
            default = "",
            type = "character",
            help = "Input target species name (e.g. mouse)"
        ),
        make_option(
            "--biomart_url",
            default = "",
            type = "character",
            help = "Biomart URL"
        ),
        make_option(
            "--release_version",
            default = "",
            type = "character",
            help = "Ensembl version number"
        ),
        make_option(
            "--output",
            default = "",
            type = "character",
            help = "Output filename"
        )
    )
    opt <- parse_args(OptionParser(option_list = option_list),
        args = commandArgs(trailingOnly = TRUE)
    )
    return(opt)
}
args <- get_args()
source_genes_biomart <- args$source_genes
target_genes_biomart <- args$target_genes
source_latin_prefix <- args$source_latin_prefix
target_latin_prefix <- args$target_latin_prefix
source_name <- args$source_name
target_name <- args$target_name
biomart_gene_dataset <- args$biomart_gene_dataset
biomart_url <- args$biomart_url
release_version <- args$release_version

output_filename <- args$output

target_gene_ensembl <- paste0(target_latin_prefix, "_gene_ensembl")

target_mart <- useEnsembl(
    biomart = "ENSEMBL_MART_ENSEMBL", host = biomart_url,
    dataset = target_gene_ensembl, version = release_version
)
message("Read BioMart files")
target_anno <- read.delim(target_genes_biomart)
source_anno <- read.delim(source_genes_biomart)

source_homolog_associated_gene_name <- paste0(source_latin_prefix, "_homolog_associated_gene_name")
source_homolog_orthology_type <- paste0(source_latin_prefix, "_homolog_orthology_type")

message(paste(c("Get", source, "orthologs")))
source_orthologs <- getBM(
    attributes = c(
        "external_gene_name",
        source_homolog_associated_gene_name,
        source_homolog_orthology_type
    ),
    filters = "external_gene_name", values = subset(target_anno, gene_biotype == "protein_coding")$gene_name, mart = target_mart
)

message("First filter on orthologs being ortholog_one2one")
target_source_orthologs <- subset(source_orthologs, source_orthologs[, source_homolog_orthology_type] == "ortholog_one2one",
    select = c("external_gene_name", source_homolog_associated_gene_name)
)

message("Keep rows with an external gene name and protein coding orthologs")
target_source_orthologs <- subset(
    target_source_orthologs,
    external_gene_name != "" &
        target_source_orthologs[, source_homolog_associated_gene_name] %in% subset(source_anno, gene_biotype == "protein_coding")$gene_name
)
target_source_orthologs <- target_source_orthologs[!duplicated(target_source_orthologs), ]
colnames(target_source_orthologs) <- c(target_name, source_name)
write.table(target_source_orthologs, output_filename, quote = F, sep = "\t", row.names = F, col.names = F)
