# ----------------------------------------------------------------------------
# Libraries and setup ----

library("Seurat")
library("tidyverse")
library("argparse")

# Command line arguments ----
parser <- ArgumentParser(description = "Filter a seurat object")

parser$add_argument('--filter', '-f', help = 'Name of filter function')
parser$add_argument('--seurat', '-s',
                    help = 'seurat object, saved as .rds')
parser$add_argument('--output', '-o', help = 'Directory to save output')

args <- parser$parse_args()

# ----------------------------------------------------------------------------
# Functions ----

# Filtering by TCR beta presence ----
tcr_beta <- function(seurat_obj){
    has_tcr <- ! is.na(seurat_obj@meta.data$CTaa2)
    seurat_obj <- subset(seurat_obj, cells = which(has_tcr))
    seurat_obj
}

one_tcr_beta <- function(seurat_obj){
    seurat_obj <- subset(seurat_obj, n_betas == 1)
    seurat_obj
}

# Filtering by CD8 expression and TCR beta presence ----
.cd8_tcr <- function(seurat_obj, fn){
    cd8_expressed <- seurat_obj$CD8B
    has_tcr <- ! (seurat_obj@meta.data$tcr_presence == "no_tcr")
    keep <- getFunction(fn)(cd8_expressed, has_tcr) 
    seurat_obj <- subset(seurat_obj, cells = which(keep))
    return(seurat_obj)
}


cd8_or_tcr <- function(seurat_obj){
    return(.cd8_tcr(seurat_obj, "|"))
}

cd8_and_tcr <- function(seurat_obj){
    return(.cd8_tcr(seurat_obj, "&"))
}

filter_seurat <- function(args){
    if (! dir.exists(args$output)){ dir.create(args$output) }
    
    seurat_obj <- read_rds(args$seurat)   
    seurat_obj <- getFunction(args$filter)(seurat_obj)
    seurat_out <- file.path(args$output, "filtered_seurat.rds")
    write_rds(seurat_obj, seurat_out)
    write_csv(seurat_obj[[]], gsub("rds$", "csv.gz", seurat_out)) # add .gz
}

# ----------------------------------------------------------------------------   
# Run ----

filter_seurat(args)
