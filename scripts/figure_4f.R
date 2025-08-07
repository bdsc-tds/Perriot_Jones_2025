# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("khroma")
library("tidyverse")
library("Seurat")
library("ComplexHeatmap")
library("RColorBrewer")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', '-s',
                    help = 'Seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
parser$add_argument('--workdir',  '-w', 
                    help = 'Working directory')

args <- parser$parse_args()

# ----------------------------------------------------------------------------
# Functions ----

# blue_yellow_red palette ----
blue_yellow_red <- function(){
    return(rev(brewer.pal(n = 7, name = "RdYlBu")))
}

pb_heatmap <- function(pseudo, markers, palette,  ...) {
    anno <- markers$cat_label[match(rownames(pseudo), markers$gene)]
    cats <- unique(anno) 
    
    if (length(cats) == 1){
        cluster_rows <- TRUE
        row_split <- NULL
        row_ha <- NULL
    } else {
        cluster_rows <- cluster_within_group(t(pseudo), anno)
        row_split <- length(cats)
        row_ha <- rowAnnotation(Category = anno,
                                col = list(Category = 
                                               structure(color("light")(length(cats)),
                                                         names = cats)),
                                show_legend = c(FALSE),
                                show_annotation_name = FALSE)
    }
    
    heatmap_args <- list(matrix = pseudo,
                         col = palette,
                         cluster_columns = TRUE,
                         show_column_dend = FALSE,
                         show_row_dend = FALSE, 
                         column_title_gp = gpar(fontsize = 10),
                         row_names_gp = gpar(fontsize = 7),
                         column_names_gp = gpar(fontsize = 7),
                         row_split = row_split,
                         row_title_gp = gpar(fontsize = 7.4),
                         heatmap_legend_param = list(title = "Scaled \nexpression"),
                         left_annotation = row_ha,
                         cluster_rows = cluster_rows,
                         column_names_rot = 45,
                         row_gap = unit(2, "mm"))
    
    heatmap_args <- modifyList(heatmap_args, list(...))
    
    return(do.call(Heatmap, heatmap_args))
}


# reactivity analyses for single marker set ----
pb_marker_set <- function(all_clones,
                          markers,
                          palette = blue_yellow_red(),
                          group_by = "rx_by_cnd",
                          agg_method = "pseudobulk",
                          ...){
    
    # Subset to genes of interest, aggregate expression
    clones <- subset(all_clones,
                     Sample != "Ri01_5m")
    
    pseudo_cat <- AggregateExpression(clones,
                                      group.by = group_by,
                                      return.seurat = TRUE)
    
    pseudo_cat <- LayerData(pseudo_cat, "scale.data")
    
    if (group_by == "seurat_clusters"){
        colnames(pseudo_cat) <- gsub("g", "", colnames(pseudo_cat))
    }
    
    pseudo_cat <- pseudo_cat[intersect(markers$gene, rownames(pseudo_cat)), ]
    
    h <- pb_heatmap(pseudo_cat,
                    markers,
                    palette = palette,
                    ...)
    return(h)
    
}

# main ----
main <- function(args){
    
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
    h <- pb_marker_set(all_clones = seurat_subs, 
                       markers = markers,
                       group_by = "seurat_clusters",
                       column_names_rot = 0,
                       row_title_gp = gpar(fontsize = 7.4),
                       column_title_gp = gpar(fontsize = 10),
                       column_names_gp = gpar(fontsize = 10))
    
    pdf(file.path(args$results, "fig_4f.pdf"), width = 5.3, height = 8)
    print(h)
    dev.off()
}

# ----------------------------------------------------------------------------
main(args)
