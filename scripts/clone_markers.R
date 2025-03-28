# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("janitor")
library("tidyverse")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--clones', help = 'Seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
parser$add_argument('--workdir',
                    help = "Working directory, for loading scripts")

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/markers_sp1.R"))
source(file.path(args$workdir, "scripts/funcs_custom_heatmaps.R"))

# ----------------------------------------------------------------------------
# Make pseudobulk ----
make_pseudobulk <- function(obj, idents){
    gp_by = unique(c("Sample", "beta_aa", idents))
    clone_pseudo <- AggregateExpression(obj,
                                        group.by = gp_by,
                                        return.seurat = TRUE)
    
    data_layer <- LayerData(clone_pseudo, "data")
    na_col <- apply(data_layer, 2, function(x) all(is.na(x)))
    data_layer[, names(na_col)[na_col]] <- 0
    LayerData(clone_pseudo, "data") <- data_layer
    clone_pseudo <- ScaleData(clone_pseudo, features = Features(clone_pseudo))
    Idents(clone_pseudo) <- idents
    return(clone_pseudo)
}

# Dotplot of cluster markers ----
make_dotplot <- function(seurat_obj, out_fname, ht = 5, scale = TRUE){
    
    pdf(out_fname, height = ht)
    p <- DotPlot(seurat_obj,
                 scale = scale,
                 features = Features(seurat_obj)) + 
        coord_flip() +
        scale_color_viridis_c() +
        labs(y = NULL, x = NULL) +
        theme(panel.grid = element_line(color = "gray")) +
        guides(size = guide_legend(title = "Percent\nExpressed"),
               colour = guide_legend(title = "Average\nExpression"))
    print(p)
    dev.off()
}

# main ----
main <- function(args){
    
    # Setup ----
    results <- file.path(args$results, "dotplots")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }
    
    marker_sets <- functional_markers()
    
    seurat_obj <- read_rds(args$clones)
    seurat_obj <- subset(seurat_obj, Sample != "Ri01_5m")
    
    seurat_obj[[]] <- seurat_obj[[]] %>%
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
    # Dotplots for marker lists ----
    template <- file.path(results, "%s_%s.pdf")
    
    for (mk_nm in names(marker_sets)){
        print(mk_nm)
        markers <- marker_sets[[mk_nm]]
        if (length(intersect(markers, Features(seurat_obj))) == 0){
            print(sprintf("no markers found: %s", mk_nm))
            next()
        }
        
        out_nm <- gsub("[-\\/ ]", "_", mk_nm)
        marker_subset <- subset(seurat_obj, features = markers)
        DefaultAssay(marker_subset) <- "RNA"
        marker_subset <- ScaleData(marker_subset, features = markers)
        Idents(marker_subset) <- "rx_by_cnd"
    
        make_dotplot(marker_subset, 
                     sprintf(template, "dotplot_unscaled_no_5m", out_nm),
                     scale = FALSE)
        make_dotplot(marker_subset, 
                     sprintf(template, "dotplot_scaled_no_5m", out_nm),
                     scale = TRUE)
    }
    
    # Average expression for reactivity x condition ----
    tested <- subset(seurat_obj, reactive == TRUE | reactive == FALSE)
    
    Idents(tested) <- "rx_by_cnd"
    avg <- AverageExpression(tested, assays = "RNA")[["RNA"]] %>%
        as_tibble(rownames = "gene")
    rx_res <- file.path(args$results,
                        "tables/reactivity_by_condition_avg_expr.csv")
    write_csv(avg, rx_res)
    
    # Average expression for clones ----
    Idents(tested) <- "tcr_name"
    avg <- AverageExpression(tested, assays = "RNA")[["RNA"]] %>%
        as_tibble(rownames = "gene") 
    colnames(avg) <- gsub("-", "_", colnames(avg))
    clone_res <- file.path(args$results,
                        "tables/tested_clones_avg_expr.csv")
    write_csv(avg, clone_res)
    
    # UMAP for clones ----
    umap_res <- file.path(args$results, "umap")

    if (! file.exists(umap_res)) { dir.create(umap_res, recursive = TRUE) }
    
    tested <- FindVariableFeatures(tested)
    tested <- ScaleData(tested)
    tested <- RunPCA(tested)
    
    pdf(file.path(umap_res, "pca_dim_red.pdf"), width = 10, height = 20)
    p <- DimHeatmap(tested, dims = 1:15, balanced = TRUE)
    print(p)
    dev.off()
    
    pdf(file.path(umap_res, "elbow.pdf"))
    p <- ElbowPlot(tested)
    print(p)
    dev.off()
    
    tested <- RunUMAP(tested, assay = "RNA", dims = 1:15)
    
    pdf(file.path(umap_res, "umap.pdf"))
    p <- DimPlot(tested, group.by = "rx_by_cnd") +
        labs(x = "UMAP 1", y = "UMAP 2", title = NULL)
    print(p)
    dev.off()
    
    #de_genes <- read_csv()
    
    
    # Average expression and percentage of cells expressing list 1_SP ----
    sp <- markers_1_sp()
    sp_obj <- subset(seurat_obj, features = sp)

    # Get percentage and average expression from dotplot data
    Idents(sp_obj) <- "rx_by_cnd"
    sp_obj <- ScaleData(sp_obj, features = sp)
    
    dp_dat <- DotPlot(object = sp_obj,
                      features = sp,
                      group.by = "rx_by_cnd")$data %>%
        dplyr::select(`features.plot`, `avg.exp`, `pct.exp`, id) %>%
        janitor::clean_names() %>%
        tidyr::pivot_wider(id_cols = features_plot,
                           names_from = id,
                           values_from = c(avg_exp, pct_exp)) %>%
        janitor::clean_names() %>%
        dplyr::rename(gene = features_plot) %>%
        readr::write_csv(file.path(args$results, "tables/markers_1_sp_expr.csv"))
    
    # Heatmap of expression for list 1_SP ----
    ht <- 12

        results <- file.path(args$results, "sp_heatmap")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }
    
    pdf(file.path(results, "unclustered.pdf"), height = ht)
    h <- DoHeatmap(sp_obj,
                   features = Features(sp_obj),
                   group.by = "rx_by_cnd") +
        theme(axis.text.y = element_text(size = 1))
    print(h)
    dev.off()
    
    pdf(file.path(results, "clustered.pdf"), height = ht)
    h <- heatmap_w_labs(sp_obj,
                        col_group = "rx_by_cnd", 
                        row_group = FALSE,
                        show_column_names = FALSE,
                        row_names_gp = gpar(fontsize = 1))
    print(h)
    dev.off()
    
    # Pseudobulk heatmap ----
    clone_pseudo <- make_pseudobulk(sp_obj, "rx_by_cnd")
    
    pdf(file.path(results, "pseudobulk_unclustered.pdf"), height = ht) 
    h <- DoHeatmap(clone_pseudo,
                   features = Features(sp_obj),
                   group.by = "rx_by_cnd") +
        theme(axis.text.y = element_text(size = 1))
    print(h)
    dev.off()
    
    pdf(file.path(results, "pseudobulk_clustered.pdf"), height = ht)
    h <- heatmap_w_labs(clone_pseudo,
                        col_group = "rx_by_cnd", 
                        row_group = FALSE,
                        row_names_gp = gpar(fontsize = 1),
                        column_names_gp = gpar(fontsize = 6))
    print(h)
    dev.off()
}


