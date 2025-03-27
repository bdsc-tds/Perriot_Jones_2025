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

# _______________________________________________________
get_markers <- function(){
    return(structure(
        c("CD74", "HSPA8", "HSP90AA1", "HLA-DRB5", "HSP90AB1", "HSPA5", 
          "TAP1", "IFI30", "TNF", "HLA-DMA", "IFNG", "HLA-DPB1", "HLA-DRA", 
          "CALR", "HLA-DQA1", "HLA-DRB1", "HLA-DPA1", "HLA-DQB1", "CDKN1A", 
          "SPI1", "ICAM1", "RELB", "ZFP36", "EGR1", "EGR2", "JUN", "MAP2K2", 
          "TGFB1", "FOS", "NFKBIA", "VDAC2", "SLC25A5", "RAN", "MAP2K3", 
          "GADD45B", "GADD45A", "TNFAIP3", "IRF7", "VIM", "CD44", "IRF9", 
          "CSF1", "CCL3", "ATP6V0C", "RARA", "IL21R")))
    
}


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
    markers <- get_markers()
    
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