# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ComplexHeatmap")
library("ggplot2")
library("tidyverse")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', '-s',
                    help = 'Seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
parser$add_argument('--clones', help = 'Clone table')
parser$add_argument('--workdir',
                    help = "Working directory, for loading scripts")

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/funcs_EnhancedVolcano_mod.R"))

# ----------------------------------------------------------------------------

# filter_tcr_genes ----
filter_tcr_genes <- function(markers){
    markers %>%
        dplyr::filter(! grepl("^TR[ABG][VDJ]", gene))
} 


# purple_and_yellow ----
purple_and_yellow <- function(disp.min = -2.5, disp.max = 2.5){
    pp_yl <- (disp.max - disp.min)/49
    pp_yl <- seq(disp.min, disp.max, by = pp_yl) 
    pp_yl_pal <- circlize::colorRamp2(pp_yl, PurpleAndYellow())
    return(pp_yl_pal)
}

# Heatmap ----
heatmap_w_labs <- function(obj,
                           col_labs = tcr_cat_to_label(),
                           disp_min = -2.5,
                           disp_max = 2.5,
                           col_group = as.factor(obj$tcr_category),
                           row_group = TRUE,
                           column_title_gp = gpar(fontsize = 10),
                           column_names_gp = gpar(fontsize = 10),
                           row_names_gp = gpar(fontsize = 3),
                           ...){
    
    pp_yl_pal <- purple_and_yellow(disp.min = disp_min, disp.max = disp_max) 
    dat <- LayerData(obj, "scale.data")
    
    col_group <- col_labs[obj[[col_group]][,1]]
    if (! isTRUE(row_group)) { 
        row_group <- cluster_within_group(t(dat), row_group)
    }
    
    Heatmap(dat,
            cluster_columns = cluster_within_group(dat, col_group),
            cluster_rows = row_group,
            show_column_dend = FALSE,
            show_row_dend = FALSE, 
            column_split = length(unique(obj$tcr_category)),
            column_title_gp = column_title_gp,
            column_names_gp = column_names_gp, 
            row_names_gp = row_names_gp,
            column_labels = gsub("_.*", "", colnames(dat)),
            col = pp_yl_pal,
            heatmap_legend_param = list(title = "Scaled \nexpression"),
            ...)
}

# Marker lists ----
get_markers <- function(){
    list("NK-like" = 
             c("KLRC1", "KLRC2", "KLRG1", "KLRK1", "KIR2DL1", "KIR2DL2", "KIR2DL3", 
               "KIR2DL4", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR3DL1", 
               "KIR3DL2", "KIR3DS1", "KIR3DS2"),
         "Activation" = 
             c("CX3CR1", "TFRC", "HLA-DRA", "HLA-DRB1", "IL2RA", "IL2RB", 
               "CD44", "ITK", "NFATC2", "NFKB1", "NME1", "CD27", "CD28",
               "ICOS", "CD40LG", "TNFRSF4", "TNFRSF9"),
         "Effector function/Cytoxicity" = 
             c("FASLG", "IFNG", "TNF", "PRF1", "GZMA", "GZMB", "GZMK", "GZMM"),
         "Exhaustion" = 
             c("PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "CD244", "CD160", 
               "TOX", "EOMES", "BATF", "TBX21", "NR4A1", "PRDM1", "ENTPD1"),
         "NK-like" = 
             c("KLRC1", "KLRC2", "KLRG1", "KLRK1", "KIR2DL1", "KIR2DL2", "KIR2DL3", 
               "KIR2DL4", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR3DL1", 
               "KIR3DL2", "KIR3DS1", "KIR3DS2"),
         "Proliferation/survival" = 
             c("MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MYBL2", "MYC", "CDK1", 
               "CDK2", "CDK4", "CDK6", "CCNA1", "CCNB1", "CCND1", "CCND2", "CCND3", 
               "CCNE1", "E2F1", "FOXM1", "BUB1", "TOP2A", "MKI67", "PCNA", "PLK1", 
               "BCL2", "BCL6", "TCF7"),
         "Tissue residency/migration" =
             c("CD69", "ITGA1", "ITGA4", "ITGAE", "CXCR3", "CXCR4", "CCR5", 
               "CXCR6", "CCR7", "S1PR1", "SELL")
    )
}

# Make single dotplot ----  
make_dotplot <- function(seurat_obj, markers, out_fname,
                         ht = 5, scale = TRUE){
    
    pdf(out_fname, height = ht)
    p <- DotPlot(seurat_obj,
                 scale = scale,
                 features = Features(seurat_obj)) + 
        coord_flip() +
        scale_color_gradient2(low="lightgray",
                              mid="#E1C7C2",
                              high="#e60000") +
        labs(y = NULL, x = NULL) +
        theme(panel.grid = element_line(color = "gray")) +
        guides(size = guide_legend(title = "Percent\nExpressed"),
               colour = guide_legend(title = "Average\nExpression"))
    print(p)
    dev.off()
}

# Make dotplots with and without scaling  ----  
make_dotplots <- function(seurat_obj, results, marker_sets){
    template <- file.path(results, "%s_%s.pdf")
    
    for (mk_nm in names(marker_sets)){
        markers <- marker_sets[[mk_nm]]
        if (length(intersect(markers, Features(seurat_obj))) == 0){
            print(sprintf("no markers found: %s", mk_nm))
            next()
        }
        
        out_nm <- gsub("[-\\/ ]", "_", mk_nm)
        marker_subset <- subset(seurat_obj, features = markers)
        DefaultAssay(marker_subset) <- "RNA"
        marker_subset <- ScaleData(marker_subset, features = markers)
        Idents(marker_subset) <- as.numeric(marker_subset$seurat_clusters)
        
        make_dotplot(marker_subset, markers,
                     sprintf(template, "dotplot_unscaled", out_nm),
                     scale = FALSE)
        make_dotplot(marker_subset, markers,
                     sprintf(template, "dotplot_scaled", out_nm),
                     scale = TRUE)
    }
}


# Filter DE table ----
filter_de <- function(de, max_n, pval_min = 0.05){
    
    n_sig <- sum(de$p_val_adj <= pval_min)
    if (n_sig >= max_n){
        de <- de %>%
            dplyr::filter(p_val_adj <= pval_min) %>%
            dplyr::arrange(desc(abs(avg_log2FC))) %>%
            dplyr::slice_head(n = max_n)
    }
    de
}

# Make heatmaps with and without tcr gene filtering ----
make_heatmaps <- function(obj, de, out_fname, max_n, remove_tcr = FALSE,
                         pval_min = 0.05, wd = 7, ht = 7, ...){
    if (isTRUE(remove_tcr)){ 
        de <- filter_tcr_genes(de)
        out_fname <- gsub("(\\.pdf|\\.png)$", "_no_tcr(\\1)", out_fname)    
    }
    
    n_sig <- sum(de$p_val_adj <= pval_min)
    if (n_sig >= max_n){
        de <- de %>%
            dplyr::filter(p_val_adj <= pval_min) %>%
            dplyr::arrange(desc(abs(avg_log2FC))) %>%
            dplyr::slice_head(n = max_n)
    }

    print(nrow(de))
    
    pdf(out_fname, width = wd, height = ht)
    h <- DoHeatmap(obj, features = de$gene, ...)
    print(h)
    dev.off()
    
    # TO DO: custom heatmaps split by samples
}

# Run differential expression analyses ----
run_diff_expr <- function(obj, results_dir, idents, ident_1, ident_2,
                          name, markers = get_markers(), max_n = 100,
                          pval_min = 0.05, wd = 7, ht = 7, ...){
    res_dir <- file.path(results_dir, name)
    if (! file.exists(res_dir)) { dir.create(res_dir, recursive = TRUE) }
    
    obj <- ScaleData(obj, features = Features(obj))
    
    # Run diff expr
    Idents(obj) <- obj[[idents]][, 1]
    de <- FindMarkers(obj, ident.1 = ident_1, ident.2 = ident_2, ...) %>%
        dplyr::as_tibble(rownames = "gene") %>%
        relocate(gene)
    
    # CSV of results
    write_csv(de, file.path(res_dir, "diff_expr_full.csv"))
    
    # Volcano of results
    pdf(file.path(res_dir, "volcano.pdf"))
    volc_labels <- filter_de(de, 20) %>% dplyr::pull(gene)
    p <- volc_mod(de,
                  xlab = expression("Average " * log[2] * "FC"),
                  ylab = expression(-log[10]*"(adjusted p-value)"),
                  labels = volc_labels,
                  labelsFrom = "gene",
                  x = "avg_log2FC",
                  y = "p_val_adj")
    print(p)
    dev.off()
    
    tryCatch({dev.off()}, error = function(e){return()})
    
    # Heatmap
    heatmap_out <- file.path(res_dir, "heatmap.pdf")
    print(heatmap_out)
    make_heatmaps(obj, de, heatmap_out, max_n,
                  pval_min, wd, ht, group.by = idents)
    
    # Customised heatmap
    #pdf(file.path(res_dir, "heatmap_clustered.pdf"))
    #h <- heatmap_w_labs(obj,
    #                    col_labs = structure(levels(Idents(obj)),
    #                                         names = levels(Idents(obj))), 
    #                    col_group = idents, 
    #                    row_group = FALSE)
    #print(h)
    #dev.off()
    
    
    # Pseudobulk
    #obj_pseudo <- AggregateExpression(obj, by = c("sample", "beta_aa"))
    #de_pseudo <- FindMarkers(obj, ident.1 = ident_1, ident.2 = ident_2,
    #                         test.use = "DESeq2")
    
    
    # Dot plots for marker set H
    make_dotplots(obj, res_dir, markers)
    
}


# Filter clones for minimum number of cells ---- 
clones_min_n <- function(md, min_cells = 5){
    keep_clones <- md %>%
        dplyr::group_by(Sample, beta_aa) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::filter(n >= min_cells) %>%
        dplyr::select(-n)
    md %>% semi_join(keep_clones)
}

# Get the cluster with the most clones of interest ----
get_cluster_coi <- function(md){
    md %>%
        dplyr::filter(is_coi == "coi") %>%
        dplyr::group_by(seurat_clusters) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::arrange(desc(n)) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::pull(seurat_clusters)
}

# main ----
main <- function(args, min_cells = 5, ...){
    
    # Setup ---- 
    
    if (! file.exists(args$results)) { 
        dir.create(args$results, recursive = TRUE)
    }
    
    seurat_obj <- read_rds(args$seurat)
    selected_clones <- read_csv(args$clones) %>% unique()
        
    # This only works as data was filtered for exactly one tcr beta
    seurat_obj[[]] <- seurat_obj[[]] %>%
        dplyr::mutate(cdr3_aa_beta = CTaa2) %>%
        tidyr::separate(TCR2,
                        into = c("trbv", "trbd", "trbj", "trbc"),
                        sep = "\\.",
                        remove = FALSE) %>%
        dplyr::left_join(selected_clones,
                         by = c("Sample","cdr3_aa_beta", "trbv", "trbj"),
                         relationship = "many-to-one")

    coi_cluster <- get_cluster_coi(seurat_obj[[]])
      
    run_de <- purrr::partial(run_diff_expr, results = args$results, ...)
    
    # Subset to clones where reactivity information is known
    clones <- subset(seurat_obj, reactive %in% c(TRUE, FALSE))
    write_rds(clones, file.path(dirname(args$seurat), "reactive_clones.rds"))
    
    clone_pseudo <- AggregateExpression(clones,
                                        group.by = c("Sample", "beta_aa"),
                                        return.seurat = TRUE)
    write_rds(clones,
              file.path(dirname(args$seurat),
                        "reactive_clones_pseudobulk.rds"))
    
    # Differential expression -----
    # USE Ri_01 dis, remove Ri_01 5m
    
    # Subset to just reactive clones, test HD v Ri
    rx <- subset(clones, reactive == "TRUE")
    run_de(rx, "condition", "Ri", "HD", "reactive_Ri_v_HD")
    
    # Subset to Ri, test reactive verus non-reactive
    # CHECK NOT Ri_01 dis??
    ri <- subset(clones, condition == "Ri")
    run_de(ri, "reactive", "TRUE", "FALSE", "Ri_rx_v_non_rx")
    
    # Subset to HD, test reactive verus non-reactive
    hd <- subset(clones, condition == "HD")
    run_de(hd, "reactive", "TRUE", "FALSE", "HD_rx_v_non_rx")
    
    # Reactive verus non-reactive, all_samples
    run_de(clones, "reactive", "TRUE", "FALSE", "all_samples_rx_v_non_rx")
    
    # Subset to cluster of interest, clones with at least min_cells
    coi_cluster <- subset(seurat_obj,
                          seurat_clusters == coi_cluster)
    keep_clones <- clones_min_n(coi_cluster[[]], min_cells)
    
    
    # To do: subset by n cells
    
    #for (mk_nm in names(marker_sets)){
    #    marker_subset <- subset(seurat_obj, features = markers)
    #    DefaultAssay(marker_subset) <- "RNA"
    #    marker_subset <- ScaleData(marker_subset, features = markers)
    #    Idents(marker_subset) <- as.numeric(marker_subset$seurat_clusters)
    #}
    
}

# ----------------------------------------------------------------------------
main(args)
