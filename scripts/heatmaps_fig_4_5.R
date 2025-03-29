# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("tidyverse")
library("scales")
library("Seurat")
library("ComplexHeatmap")
library("khroma")
library("viridis")
library("RColorBrewer")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', '-s',
                    help = 'Seurat object')
parser$add_argument('--clones', '-c',
                    help = 'Tested clones, rds')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
parser$add_argument('--workdir',  '-w', 
                    help = 'working directory')

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/funcs_custom_heatmaps.R"))

# ----------------------------------------------------------------------------
# Functions ----

# cluster_heatmap ----
cluster_heatmap <- function(pseudo_cl, markers, palette = blue_and_yellow()){
    dat <- LayerData(pseudo_cl, "scale.data")
    anno <- markers$cat_label[match(Features(pseudo_cl), markers$gene)]
    cats <- unique(anno) 
    row_ha = rowAnnotation(Category = anno,
                           col = list(Category = 
                                          structure(color("light")(length(cats)),
                                                    names = cats)),
                           show_legend = c(FALSE),
                           show_annotation_name = FALSE)
    
    h <- Heatmap(dat,
            cluster_columns = TRUE,
            show_column_dend = FALSE,
            show_row_dend = FALSE, 
            column_title_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            row_names_gp = gpar(fontsize = 7),
            row_split = length(unique(anno)),
            row_title_gp = gpar(fontsize = 7.4),
            column_labels = gsub("g", "", colnames(dat)),
            col = palette,
            heatmap_legend_param = list(title = "Scaled \nexpression"),
            left_annotation = row_ha,
            cluster_rows = cluster_within_group(t(dat), anno),
            column_names_rot = 0)
    return(h)
}

# reactivity heatmap ----
rx_heatmap <- function(pseudo_cat, markers, palette = blue_and_yellow()) {
    anno <- markers$cat_label[match(Features(pseudo_cl), markers$gene)]
    cats <- unique(anno) 
    row_ha <- rowAnnotation(Category = anno,
                           col = list(Category = 
                                          structure(color("light")(length(cats)),
                                                    names = cats)),
                           show_legend = c(FALSE),
                           show_annotation_name = FALSE)
    
    Heatmap(LayerData(pseudo_cat, "scale.data"), 
            cluster_columns = TRUE,
            show_column_dend = FALSE,
            show_row_dend = FALSE, 
            column_title_gp = gpar(fontsize = 10),
            row_names_gp = gpar(fontsize = 7),
            column_names_gp = gpar(fontsize = 7),
            row_split = length(unique(anno)),
            row_title_gp = gpar(fontsize = 7.4),
            column_labels = unique(pseudo_cat$orig.ident),
            col = palette,
            heatmap_legend_param = list(title = "Scaled \nexpression"),
            left_annotation = row_ha,
            cluster_rows = cluster_within_group(t(dat), anno),
            column_names_rot = 45,
            row_gap = unit(2, "mm"))
}

# main ----
main <- function(args){
    results <- file.path(args$results, "heatmaps_fig_4_5")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }
    bu_yl_dir <- file.path(results, "blue_and_yellow")
    if (! file.exists(bu_yl_dir)) { dir.create(bu_yl_dir) }
    viridis_dir <- file.path(results, "viridis")
    if (! file.exists(viridis_dir)) { dir.create(viridis_dir) }
    bu_yl_rd_dir <- file.path(results, "blue_yellow_red")
    if (! file.exists(bu_yl_rd_dir)) { dir.create(bu_yl_rd_dir) }
    
    bu_yl_rd_pal <- rev(brewer.pal(n = 7, name = "RdYlBu"))
    
    markers <- read_csv(file.path(args$workdir,
                                  "data/processed/gene_lists_4_5.csv")) %>%
        dplyr::mutate(cat_label = gsub(" ", "\n", category),
                      cat_label = gsub("\\/", " \\/\n", cat_label))
    
    seurat_obj <- read_rds(args$seurat)
    
    # Exclude Ri01_5m 
    seurat_subs <- subset(seurat_obj, Sample != "Ri01_5m",
                          features = markers$gene)
    seurat_subs <- ScaleData(seurat_subs, features = markers$gene) 
    Idents(seurat_subs) <- "seurat_subs"
        
    write_rds(seurat_subs,
              file.path(dirname(args$seurat), "fig_4_5_subset.rds"))

    seurat_subs <- read_rds("data/processed/one_tcr_beta/fig_4_5_subset.rds")
    
    # Aggregate across clusters ----
    pseudo_cl <- AggregateExpression(seurat_subs,
                                     group.by = "seurat_clusters",
                                     return.seurat = TRUE)
    
    pdf(file.path(bu_yl_dir, "category_by_cluster.pdf"),
        width = 5.3, height = 8)
    h <- cluster_heatmap(pseudo_cl, markers)
    print(h)
    dev.off()
    
    pdf(file.path(viridis_dir, "category_by_cluster.pdf"),
        width = 5.3, height = 8)
    h <- cluster_heatmap(pseudo_cl, markers, palette = viridis(100))
    print(h)
    dev.off()
    
    pdf(file.path(bu_yl_rd_dir, "category_by_cluster.pdf"),
        width = 5.3, height = 8)
    h <- cluster_heatmap(pseudo_cl, markers, palette = bu_yl_rd_pal)
    print(h)
    dev.off()
    
    
    rm(seurat_subs)
    
    # Aggregate across reactivity categories 
    
    #clones <- read_rds("data/processed/one_tcr_beta/reactive_clones.rds")
    clones <- read_rds(args$clones)
    clones <- subset(clones, Sample != "Ri01_5m", features = markers$gene)
    
    clones[[]] <- clones[[]] %>%
        dplyr::mutate(condition = case_when(condition == "Ri01_dis" ~ "Ri",
                                            TRUE ~ condition),
                      rx_by_cnd =
                          case_when(condition == "HD" & reactive == "TRUE" ~
                                        "HD reactive",
                                    condition == "HD" & reactive == "FALSE" ~
                                        "HD non-reactive",
                                    condition == "Ri" & reactive == "TRUE" ~
                                        "Ri reactive",
                                    condition == "Ri" & reactive == "FALSE" ~
                                        "Ri non-reactive"))
    
    pseudo_cat <- AggregateExpression(clones,
                                      group.by = "rx_by_cnd",
                                      return.seurat = TRUE)
    
    pdf(file.path(bu_yl_dir, "category_by_reactivity.pdf"),
        width = 3, height = 8)
    h <- rx_heatmap(pseudo_cat, markers)
    print(h)
    dev.off()
    
    pdf(file.path(viridis_dir, "category_by_reactivity.pdf"),
        width = 3, height = 8)
    h <- rx_heatmap(pseudo_cat, markers, palette = viridis(100))
    print(h)
    dev.off()
    
    pdf(file.path(bu_yl_rd_dir, "category_by_reactivity.pdf"),
        width = 3, height = 8)
    h <- rx_heatmap(pseudo_cat, markers, palette = bu_yl_rd_pal)
    print(h)
    dev.off()
    
}

# ----------------------------------------------------------------------------
main(args)
