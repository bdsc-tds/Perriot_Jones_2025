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

source(file.path(args$workdir, "scripts/funcs_custom_heatmaps.R"))
source(file.path(args$workdir, "scripts/markers_sp1.R"))

# _______________________________________________________

clones_min_n <- function(md, min_cells = 5){
    keep_clones <- md %>%
        dplyr::group_by(Sample, beta_aa) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::filter(n >= min_cells) %>%
        dplyr::select(-n)
    md %>% semi_join(keep_clones)
}


main <- function(args){
    results <- file.path(args$results, "string")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }
    
    wd <- 10
    ht <- 10
    markers <- string_markers()
    
    seurat_obj <- read_rds(args$seurat)
    seurat_obj <- subset(seurat_obj, features = markers)
    Idents(seurat_obj) <- "seurat_clusters"
    seurat_obj <- ScaleData(seurat_obj, features = markers)
    
    pdf(file.path(results, "heatmap_unclustered.pdf"))
    h <- DoHeatmap(seurat_obj, features = markers,
                   size = 4, angle = 0)
    print(h)
    dev.off()
    
    column_ha = HeatmapAnnotation(Sample = seurat_obj$Sample)
    
    pdf(file.path(results, "heatmap_clustered.pdf"), height = ht, width = wd)
    h <- heatmap_w_labs(seurat_obj,
                        col_group = "seurat_clusters", 
                        row_group = FALSE,
                        show_column_names = FALSE,
                        top_annotation = column_ha) # Add sample key 
    print(h)
    dev.off()
    
    keep_clones <- clones_min_n(seurat_obj, args$min_cells) ####
    seurat_subs <- subset(seurat_obj,
                          cells = rownames(keep_clones))
    
    pseudo <- AggregateExpression(seurat_subs,
                                  group.by = c("Sample",
                                               "beta_aa",
                                               "seurat_clusters"),
                                  return.seurat = TRUE)
    Idents(pseudo) <- "seurat_clusters"
    
    pdf(file.path(results, "heatmap_unclustered_pseudobulk.pdf"))
    h <- DoHeatmap(pseudo, features = markers)
    print(h)
    dev.off()
    
    pseudo_ha = HeatmapAnnotation(Sample = seurat_obj$Sample)
    
    pdf(file.path(results, "pseudobulk_clustered.pdf"), height = ht)
    h <- heatmap_w_labs(pseudo,
                        col_group = "seurat_clusters", 
                        row_group = FALSE,
                        show_column_names = FALSE,
                        top_annotation = pseudo_ha) # Add sample key 
    print(h)
    dev.off()
}

# ----------------------------------------------------------------------------
main(args)