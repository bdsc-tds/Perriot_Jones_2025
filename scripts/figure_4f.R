# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("tidyverse")
library("scales")
library("Seurat")
library("ComplexHeatmap")
library("khroma")
library("RColorBrewer")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', '-s',
                    help = 'Seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
parser$add_argument('--workdir',  '-w', 
                    help = 'working directory')

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/funcs_custom_heatmaps.R"))

# ----------------------------------------------------------------------------
# Functions ----

# main ----
main <- function(args){
    
    # Setup with all cells, genes of interest ----
    
    if (! file.exists(args$results)) { dir.create(args$results, recursive=TRUE) }

    markers <- read_csv(file.path(args$workdir,
                                  "data/processed/gene_lists_4_5.csv")) %>%
        dplyr::mutate(cat_label = gsub(" ", "\n", category),
                      cat_label = gsub("\\/", " \\/\n", cat_label))
    
    fig_4_5_subset <- file.path(dirname(args$seurat), "fig_4_5_subset.rds") 
    
    if (file.exists(fig_4_5_subset)){ # If rerunning, load saved object
        seurat_subs <- read_rds(fig_4_5_subset)
    } else { # Otherwise, generate it
        seurat_obj <- read_rds(args$seurat)
        
        # Exclude Ri01_5m remission sample
        seurat_subs <- subset(seurat_obj, Sample != "Ri01_5m",
                              features = markers$gene)
        seurat_subs <- ScaleData(seurat_subs, features = markers$gene) 
        Idents(seurat_subs) <- "seurat_subs"
        write_rds(seurat_subs, fig_4_5_subset)
    }
    
           
    # Make heatmap of expression aggregated across clusters ----
    pb_marker_partial <- purrr::partial(pb_marker_set,
                                        all_clones = seurat_subs, 
                                        markers = markers,
                                        group_by = "seurat_clusters",
                                        width = 5.3,
                                        height = 8,
                                        column_names_rot = 0,
                                        row_title_gp = gpar(fontsize = 7.4),
                                        column_title_gp = gpar(fontsize = 10),
                                        column_names_gp = gpar(fontsize = 10))
    
    pb_marker_partial(name = "category_by_cluster.pdf") # Sums
}

# ----------------------------------------------------------------------------
main(args)
