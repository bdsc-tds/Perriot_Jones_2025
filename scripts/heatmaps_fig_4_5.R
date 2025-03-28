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
    #DoHeatmap(pseudo_cl, group.by = "seurat_clusters")
    
    dat <- LayerData(pseudo_cl, "scale.data")
    #pp_yl_pal <- purple_and_yellow(disp.min = -2.5, disp.max = 2.5) 
    bu_yl_pal <- blue_and_yellow()
    
    anno <- markers$cat_label[match(Features(pseudo_cl), markers$gene)]
    cats <- unique(anno) 
    row_ha = rowAnnotation(Category = anno,
                           col = list(Category = 
                                          structure(color("light")(length(cats)),
                                          names = cats)),
                           show_legend = c(FALSE),
                           show_annotation_name = FALSE)
    
    pdf(file.path(results, "category_by_cluster.pdf"),
        width = 5.3, height = 8)
    Heatmap(dat,
            cluster_columns = TRUE,
            show_column_dend = FALSE,
            show_row_dend = FALSE, 
            column_title_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            row_names_gp = gpar(fontsize = 7),
            row_split = length(unique(anno)),
            row_title_gp = gpar(fontsize = 7.4),
            column_labels = gsub("g", "", colnames(dat)),
            col = bu_yl_pal, #pp_yl_pal,
            heatmap_legend_param = list(title = "Scaled \nexpression"),
            left_annotation = row_ha,
            cluster_rows = cluster_within_group(t(dat), anno),
            column_names_rot = 0)
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
    
    pdf(file.path(results, "category_by_reactivity.pdf"), width = 3, height = 8)
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
            col = bu_yl_pal, #pp_yl_pal,
            heatmap_legend_param = list(title = "Scaled \nexpression"),
            left_annotation = row_ha,
            cluster_rows = cluster_within_group(t(dat), anno),
            column_names_rot = 45,
            row_gap = unit(2, "mm"))
    dev.off()
}

# ----------------------------------------------------------------------------
main(args)
