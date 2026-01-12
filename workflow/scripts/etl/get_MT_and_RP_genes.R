library(data.table)
library(GenomicFeatures)
library(optparse)

###############################################################################
#' Extract MT and RP genes
#'
#' @name getSurfaceGeneList.R
#' @author Minkyung Sung, Gautier Koscielny
#' @date 7 November 2022
#' @description This script takes the surfaceome file from ethz and filter the
#' table on the `surface` keyword.
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
            "--biomart",
            default = "",
            type = "character",
            help = "Input biomart filename"
        ),
        make_option(
            "--gtf",
            default = "",
            type = "character",
            help = "Input gtf filename"
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
gtf_filename <- args$gtf
biomart_filename <- args$biomart
output_filename <- args$output

# Read the BioMart file
anno <- fread(biomart_filename)
anno <- subset(anno, gene_biotype %in% c("protein_coding", "lincRNA", "antisense"), select = c("ensembl_gene_id", "gene_name", "description"))
anno <- anno[!duplicated(anno), ]
message("makeTxDbFromGFF")
gtf_db <- makeTxDbFromGFF(gtf_filename)
message("Get genes")
# Get genes
gene_list <- as.data.frame(genes(gtf_db))
colnames(gene_list) <- c("chrom", "start", "end", "width", "strand", "ensembl_gene_id")
message("Merge Annotations...")
anno_merged <- merge(gene_list, anno, by = "ensembl_gene_id", all.x = T)
message("Define RP filters...")
RP_filt1 <- with(anno_merged, grepl("^R[Pp][SsLl]", gene_name))
RP_filt2 <- with(anno_merged, !is.na(description) & grepl("ribosomal protein", description))
message("Define MT filters...")
out <- subset(anno_merged, chrom == "MT" | RP_filt1 | RP_filt2)
message("Write output table...")
write.table(out, output_filename, quote = F, row.names = F, col.names = T, sep = "\t")
