# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("tidyverse")
library("scales")
library("Seurat")
library("ComplexHeatmap")
library("khroma")

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

# main ----
main <- function(args){
    results <- file.path(args$results, "heatmaps_fig_4_5")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }

    markers <- read_csv(file.path(args$workdir,
                                  "data/processed/gene_lists_4_5.csv"))
    
    seurat_obj <- read_rds(args$seurat)
    
    # Exclude Ri01_5m 
    seurat_subs <- subset(seurat_obj, Sample != "Ri01_5m",
                          features = markers$gene)
    seurat_subs <- ScaleData(seurat_subs, features = markers$gene) 
    Idents(seurat_subs) <- "seurat_subs"
        
    write_rds(seurat_subs,
              file.path(dirname(args$seurat), "fig_4_5_subset.rds"))

    # Aggregate across clusters ----
    pseudo_cl <- AggregateExpression(seurat_subs,
                                     group.by = "seurat_clusters",
                                     return.seurat = TRUE)
    #DoHeatmap(pseudo_cl, group.by = "seurat_clusters")
    
    dat <- LayerData(pseudo_cl, "scale.data")
    pp_yl_pal <- purple_and_yellow(disp.min = -2.5, disp.max = 2.5) 
    anno <- markers$category[match(Features(pseudo_cl), markers$gene)]
    row_ha = rowAnnotation(Category = anno,
                           col = list(Category = 
                                          structure(color("light")(7),
                                          names = unique(anno))))
    
    pdf(file.path(results, "category_by_cluster.pdf"))
    Heatmap(dat,
            cluster_columns = TRUE,
            show_column_dend = FALSE,
            show_row_dend = FALSE, 
            column_title_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            row_names_gp = gpar(fontsize = 6),
            column_labels = gsub("g", "", colnames(dat)),
            col = pp_yl_pal,
            heatmap_legend_param = list(title = "Scaled \nexpression"),
            left_annotation = row_ha,
            cluster_rows = cluster_within_group(t(dat), anno))
    dev.off()
    
    rm(seurat_subs)
    
    # Aggregate across reactivity categories 
    
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
    
    pdf(file.path(results, "category_by_reactivity.pdf"), width = 5)
    heatmap_w_labs(pseudo_cat,
                   col_group = "rx_by_cnd",
                   show_column_names = FALSE,
                   row_group = anno,
                   left_annotation = row_ha,
                   column_title_gp = gpar(fontsize = 7),
                   row_names_gp = gpar(fontsize = 6))
    dev.off()
    
}

# ----------------------------------------------------------------------------
main(args)
