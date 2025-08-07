
library("argparse")
library("ggplot2")
library("patchwork")
library("tidyverse")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', '-s',
                    help = 'Seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
args <- parser$parse_args()

# ----------------------------------------------------------------------------

main <- function(args){
    if (! file.exists(args$results)) { dir.create(args$results, recursive = TRUE) }
    
    seurat_obj <- read_rds(args$seurat)
    
    # Cluster DE, per cell
    Idents(seurat_obj) <- seurat_obj$seurat_clusters
    de <- FindAllMarkers(seurat_obj)
    write_csv(de,
              file.path(args$results, "supp_table_6_de_between_clusters.csv"))
    
}
# ----------------------------------------------------------------------------
main(args)
