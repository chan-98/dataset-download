#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)
library(readxl)

###############################################################################
#' Extract the cell surface gene list from ethz
#'
#' @name getSurfaceGeneList.R
#' @author Cedric Mendoza, Gautier Koscielny
#' @date 2 November 2022
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
            "--surfaceome",
            default = "",
            type = "character",
            help = "Input surfaceome filename"
        ),
        make_option(
            "--biomart",
            default = "",
            type = "character",
            help = "Input biomart filename"
        ),
        make_option(
            "--fantom",
            default = "",
            type = "character",
            help = "Input fantom filename"
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
surfaceome_filename <- args$surfaceome
biomart_filename <- args$biomart
fantom_filename <- args$fantom
output_filename <- args$output

# Read FANTOM matrix
# Column headings in the input file don't match the columns
# Read the header first so we can insert "hgnc_symbol" as the first column name

count_file <- file(fantom_filename, "r")
first_lines <- readLines(count_file, n = 2)
close(count_file)

header <- unlist(strsplit(first_lines[1], "\t"))
tpms <- read_tsv(fantom_filename, col_names = FALSE, skip = 1)
colnames(tpms) <- c("hgnc_symbol", header)

# Read annotation files
surface_df <- read_excel(surfaceome_filename, sheet = 2, skip = 1)
biomart_df <- read_tsv(biomart_filename)

bm_surface <- inner_join(biomart_df, surface_df, by = c("ensembl_gene_id" = "Ensembl gene"))
ff_bm_surface <- inner_join(tpms[, "hgnc_symbol"], bm_surface, by = "hgnc_symbol")

write_tsv(ff_bm_surface, file = output_filename, na = "NA", quote = "none")
# write.table(ff_bm_surface, file = output_filename, quote = FALSE, sep = "\t", col.names = TRUE)
