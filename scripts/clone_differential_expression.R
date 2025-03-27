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
parser$add_argument('--min_cells',  '-m', default = 5, 
                    help = 'Minimum cells per clone for pseudobulk')

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/funcs_EnhancedVolcano_mod.R"))

# ----------------------------------------------------------------------------

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

# filter_tcr_genes ----
filter_tcr_genes <- function(markers){
    markers %>%
        dplyr::filter(! grepl("^TR[ABGD][VDJ]", gene))
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
                           col_group,
                           col_labs = tcr_cat_to_label(),
                           disp_min = -2.5,
                           disp_max = 2.5,
                           row_group = FALSE,
                           column_title_gp = gpar(fontsize = 10),
                           column_names_gp = gpar(fontsize = 10),
                           row_names_gp = gpar(fontsize = 3),
                           show_column_names = FALSE,
                           ...){
    pp_yl_pal <- purple_and_yellow(disp.min = disp_min, disp.max = disp_max) 
    dat <- LayerData(obj, "scale.data")
    dat[, apply(dat, 2, function(x) all(is.na(x)))] <- 0
    print(range(dat))
    
    col_group <- col_labs[as.character(obj[[col_group]][,1])]
    
    if (! isTRUE(row_group) & ! isFALSE(row_group)) { 
        row_group <- cluster_within_group(t(dat), row_group)
    }
    
    Heatmap(dat,
            cluster_columns = cluster_within_group(dat, col_group),
            cluster_rows = row_group,
            show_column_dend = FALSE,
            show_row_dend = FALSE, 
            # Causing problems for one_v_all
            column_split = length(unique(col_group)), 
            column_title_gp = column_title_gp,
            column_names_gp = column_names_gp, 
            row_names_gp = row_names_gp,
            column_labels = gsub("_.*", "", colnames(dat)),
            col = pp_yl_pal,
            show_column_names = show_column_names,
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
        scale_color_viridis_c() +
        #scale_color_gradient2(low="lightgray",
        #                      mid="#E1C7C2",
        #                      high="#e60000") +
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
        #Idents(marker_subset) <- as.numeric(marker_subset$seurat_clusters)
        
        make_dotplot(marker_subset, markers,
                     sprintf(template, "dotplot_unscaled", out_nm),
                     scale = FALSE)
        make_dotplot(marker_subset, markers,
                     sprintf(template, "dotplot_scaled", out_nm),
                     scale = TRUE)
    }
}

# Filter DE table ----
# This is used for selecting genes to highlight on volcano plots
filter_de <- function(de, max_n, pval_min = 0.05){
    n_sig <- sum(de$p_val_adj <= pval_min)
    de <- de %>%
        dplyr::filter(p_val_adj <= pval_min) %>%
        dplyr::arrange(desc(abs(avg_log2FC))) %>%
        dplyr::slice_head(n = min(max_n, n_sig))
    de
}

# write top genes ----
write_top_genes <- function(de, obj, res_dir, pval_min, fc_cutoff = 0.25){
    de_genes <- de %>%
        dplyr::filter(p_val_adj <= pval_min, abs(avg_log2FC) >= fc_cutoff) %>%
        dplyr::pull(gene)
    
    strs <- gsub("\\.", "_", as.character(c(pval_min, fc_cutoff)))
    out_fname <- sprintf("de_genes_logcounts_pval_%s_fc_%s.csv",
                         strs[1], strs[2])
    
    # Log counts
    LayerData(obj,
              assay = "RNA",
              layer = "data",
              features = de_genes) %>%
        tibble::as_tibble(rownames = "gene") %>%
        write_csv(file.path(res_dir, out_fname))
    
    # Counts 
    LayerData(obj,
              assay = "RNA",
              layer = "counts",
              features = de_genes) %>%
        tibble::as_tibble(rownames = "gene") %>%
        write_csv(file.path(res_dir, gsub("logcounts", "counts", out_fname)))
    
    # Pseudobulk, no FC cutoff.
    de_pval <- de %>%
        dplyr::filter(p_val_adj <= pval_min) %>%
        dplyr::arrange(p_val_adj)
    AggregateExpression(obj, group.by = "Sample",
                        features = de_pval$gene,
                        assays = "RNA")[[1]] %>%
        tibble::as_tibble(rownames = "gene") %>%
        write_csv(file.path(res_dir, 
                            sprintf("pseudobulk_sample_pval_%s.csv", strs[1])))
}

# make volcano ----
make_volcano <- function(de, res_dir, n_labs = 20){
    print(filter_de(de, n_labs))
    volc_labels <- filter_de(de, n_labs) %>% dplyr::pull(gene)
    print("n_volc_labels")
    print(length(volc_labels))
    pdf(file.path(res_dir, "volcano.pdf"))
    p <- volc_mod(de,
                  xlab = expression("Average " * log[2] * "FC"),
                  ylab = expression(-log[10]*"(adjusted p-value)"),
                  labels = volc_labels,
                  labelsFrom = "gene",
                  x = "avg_log2FC",
                  y = "p_val_adj")
    print(p)
    dev.off()
}

# Make heatmaps ----
make_heatmaps <- function(obj, de, idents, out_fname, max_n, remove_tcr = FALSE,
                          pseudobulk = TRUE, pval_min = 0.05, wd = 7, ht = 7, ...){
    if (isTRUE(remove_tcr)){ 
        de <- filter_tcr_genes(de)
        out_fname <- gsub("(\\.pdf|\\.png)$", "_no_tcr\\1", out_fname)    
    }

    n_sig <- sum(de$p_val_adj <= pval_min)
    # If there aren't enough differentially expressed genes, ignore p_val cutoff
    if (n_sig > max_n){
        de <- de %>%
            dplyr::filter(p_val_adj <= pval_min)
    }
    
    if (nrow(de) >= max_n){
         de <- de%>%
            dplyr::mutate(pval_sig = p_val_adj <= pval_min,
                          pval_sig = factor(pval_sig, levels = c("TRUE", "FALSE")) ) %>%
            # Keep the significant genes if there are any
            dplyr::arrange(pval_sig, desc(abs(avg_log2FC)))
         print(head(de))
         de <- de %>%
            dplyr::slice_head(n = max_n)
    }
    
    obj_subs <- subset(obj, features = de$gene)
    obj_subs <- ScaleData(obj_subs, features = Features(obj_subs))
    
    # Default Seurat heatmap ----
    pdf(out_fname, width = wd, height = ht)
    h <- DoHeatmap(obj_subs, features = de$gene, ...) +
        theme(axis.text.y = element_text(size = 4))
    print(h)
    dev.off()
    
    # Customised heatmap ----

    pdf(gsub(".pdf", "_clustered.pdf", out_fname), width = wd, height = ht)
    h <- heatmap_w_labs(obj_subs,
                        col_labs = structure(levels(Idents(obj_subs)),
                                             names = levels(Idents(obj_subs))), 
                        col_group = idents, 
                        row_group = FALSE)
    print(h)
    dev.off()
    
    if (isFALSE(pseudobulk)) return()
    
    # Pseudobulk heatmap ----
    clone_pseudo <- make_pseudobulk(obj_subs, idents)
    
    pdf(gsub(".pdf", "_pseudobulk.pdf", out_fname), width = wd, height = ht)
    h <- DoHeatmap(clone_pseudo, features = de$gene, ...) +
        theme(axis.text.y = element_text(size = 4))
    print(h)
    dev.off()
    
    pdf(gsub(".pdf", "_pseudobulk_clustered.pdf", out_fname),
        width = wd, height = ht)
    h <- heatmap_w_labs(clone_pseudo,
                        col_labs = structure(levels(Idents(clone_pseudo)),
                                             names = levels(Idents(clone_pseudo))), 
                        col_group = idents, 
                        row_group = FALSE)
    print(h)
    dev.off()
    
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
        dplyr::mutate("comparison" = name) %>%
        relocate(gene, comparison)
    
    # CSV of results
    write_csv(de, file.path(res_dir, "diff_expr_full.csv"))

    de <- filter_tcr_genes(de)
    write_csv(de, file.path(res_dir, "diff_expr_rm_tcr.csv"))
    
    # Volcano of results
    make_volcano(de, res_dir, n_labs = 20)
    
    ###########
    # Stopping, just regenerating volcano
    return() #############
        
    # CSV of gene expression
    write_top_genes(de, obj, res_dir, pval_min, fc_cutoff = 0.1)

    tryCatch({dev.off()}, error = function(e){return()})
    
    # Heatmap
    heatmap_out <- file.path(res_dir, "heatmap.pdf")

    make_heatmaps(obj = obj, de = de, idents = idents, out_fname = heatmap_out,
                  max_n = max_n, pval_min = pval_min, wd = wd, ht = ht,
                  group.by = idents)
    
    # Dot plots for marker set H
    make_dotplots(obj, res_dir, markers)
    
    return(de)
}

# run_diff_expr_pb ----
run_diff_expr_pb <- function(obj, results_dir, idents, ident_1, ident_2,
                             name, markers = get_markers(), max_n = 100,
                             pval_min = 0.05, wd = 7, ht = 7, ...){
    res_dir <- file.path(results_dir, name)
    if (! file.exists(res_dir)) { dir.create(res_dir, recursive = TRUE) }
    
    pseudo <- make_pseudobulk(obj, idents)
    
    # Run diff expr
    Idents(pseudo) <- pseudo[[idents]][, 1]
    de <- FindMarkers(pseudo, ident.1 = ident_1, ident.2 = ident_2, 
                      test.use = "DESeq2", ...) %>%
        dplyr::as_tibble(rownames = "gene") %>%
        dplyr::mutate("comparison" = name) %>%
        relocate(gene, comparison)
    
    # CSV of results
    write_csv(de, file.path(res_dir, "diff_expr_pseudo_full.csv"))
    
    de <- filter_tcr_genes(de)
    write_csv(de, file.path(res_dir, "diff_expr_pseudo_rm_tcr.csv"))
    
    make_heatmaps(pseudo, de, idents,
                  out_fname = file.path(res_dir("heatmap_pseudobulk_de.pdf")),
                  max_n,
                  pseudobulk = FALSE)
    
    return(de)
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

# Make counts table ----
cluster_counts_table <- function(md, out_fname){
    md %>%
        dplyr::group_by(Sample, beta_aa, tcr_name) %>%
        dplyr::summarise(n_clone = n()) %>%
        dplyr::arrange(Sample, desc(n_clone)) %>%
        write_csv(out_fname)
}

# clones counts table ----
clones_counts_table <- function(md, out_fname){
    md %>%
        dplyr::group_by(Sample, beta_aa, tcr_name, seurat_clusters) %>%
        dplyr::summarise(n_clone = n()) %>%
        dplyr::arrange(Sample, tcr_name, desc(n_clone)) %>%
        write_csv(out_fname)
}

# main ----
main <- function(args, min_cells = 5, ...){
    
    # Setup ---- 
    if (! file.exists(args$results)) { 
        dir.create(args$results, recursive = TRUE)
        dir.create(file.path(args$results, "condition_by_reactivity"))
        dir.create(file.path(args$results, "tables"))
    }
    
    seurat_obj <- read_rds(args$seurat)
    selected_clones <- read_csv(args$clones) %>% unique()
        
    # This works as data was filtered for exactly one tcr beta
    seurat_obj[[]] <- seurat_obj[[]] %>%
        dplyr::mutate(cdr3_aa_beta = CTaa2) %>%
        tidyr::separate(TCR2,
                        into = c("trbv", "trbd", "trbj", "trbc"),
                        sep = "\\.",
                        remove = FALSE) %>%
        dplyr::left_join(selected_clones,
                         by = c("Sample","cdr3_aa_beta", "trbv", "trbj"),
                         relationship = "many-to-one")

    coi_cluster_id <- get_cluster_coi(seurat_obj[[]])
    
    # Write table of counts per cluster 
    cluster_counts_table(seurat_obj[[]] %>%
                             dplyr::filter(seurat_clusters == coi_cluster_id),
                      file.path(args$results,
                                sprintf("tabels/clone_counts_cl_%s.csv",
                                        coi_cluster_id)))
    
    run_de <- purrr::partial(run_diff_expr, results = args$results, ...)
    run_pseudo <- purrr::partial(run_diff_expr_pb, results = args$results, ...)
    
    # Subset to clones where reactivity information is known
    clones <- subset(seurat_obj, reactive %in% c(TRUE, FALSE))
    write_rds(clones, file.path(dirname(args$seurat), "reactive_clones.rds"))
    
    # Write table of reactive counts
    clones_counts_table(clones[[]],
                        file.path(args$results,
                                  "tables/reactive_clone_counts.csv"))
    
    # Differential expression -----
    clones_wo_5m <- subset(clones, Sample != "Ri01_5m")
    
    # Fix condition
    clones_wo_5m[[]] <- clones_wo_5m[[]] %>%
        dplyr::mutate(condition = case_when(condition == "Ri01_dis" ~ "Ri",
                                            TRUE ~ condition),
                      reactive_lab = case_when(reactive == "TRUE" ~ "Reactive",
                                               reactive == "FALSE" ~ "Non-reactive"))

    # Subset to just reactive clones, test HD v Ri ----
    rx <- subset(clones_wo_5m, reactive == "TRUE")
    #cnd_de <- run_de(rx, "condition", "Ri", "HD", "reactive_Ri_v_HD")
    cnd_de_pb <- run_pseudo(rx, "condition", "Ri", "HD", "reactive_Ri_v_HD_pseudo")
    
    # Subset to Ri, test reactive versus non-reactive ----
    ri <- subset(clones_wo_5m, condition == "Ri")
    #ri_de <- run_de(ri, "reactive_lab",
    #                "Reactive", "Non-reactive", "Ri_rx_v_non_rx")
    ri_de_pb <- run_pseudo(ri, "reactive_lab",
                           "Reactive", "Non-reactive", "Ri_rx_v_non_rx_pseudo")
    
    # Subset to HD, test reactive versus non-reactive ----
    hd <- subset(clones_wo_5m, condition == "HD")
    #hd_de <- run_de(hd, "reactive_lab",
    #                "Reactive", "Non-reactive", "HD_rx_v_non_rx")
    hd_de_pb <- run_pseudo(hd, "reactive_lab",
                          "Reactive", "Non-reactive", "HD_rx_v_non_rx_pseudo")
    
    # Reactive verus non-reactive, all_samples ----
    #rx_dr <- run_de(clones_wo_5m, "reactive_lab",
    #                "Reactive", "Non-reactive",
    #                "all_samples_rx_v_non_rx")
    rx_dr_pb <- run_pseudo(clones_wo_5m, "reactive_lab",
                           "Reactive", "Non-reactive",
                           "all_samples_rx_v_non_rx")
    
    # Reactive by condition ----
    clones_wo_5m[[]] <- clones_wo_5m[[]] %>%
        dplyr::mutate(rx_by_cnd = paste(condition, reactive, sep = "_"))
    Idents(clones_wo_5m) <- "rx_by_cnd"
    
    one_v_all_de <- FindAllMarkers(clones_wo_5m)
    write_csv(one_v_all_de,
              file.path(args$results, "condition_by_reactivity/de_full.csv"))
    make_heatmaps(clones_wo_5m, one_v_all_de, "rx_by_cnd",
                  max_n = 100, remove_tcr = TRUE,
                  out_fname = 
                      file.path(args$results,
                                "condition_by_reactivity/heatmap.pdf"))
    
    clones_wo_5m_pb <- make_pseudobulk(clones_wo_5m, "rx_by_cnd")
    one_v_all_de_pb <- FindAllMarkers(clones_wo_5m_pb, test.use = "DESeq2")
    write_csv(one_v_all_de_pb,
              file.path(args$results,
                        "condition_by_reactivity/de_full_pseudobulk.csv"))
    
        
    # Subset to cluster of interest, clones with at least min_cells ----
    coi_cluster <- subset(clones_wo_5m,
                          seurat_clusters == coi_cluster_id)
    keep_clones <- clones_min_n(coi_cluster[[]], args$min_cells) ####
    coi_cluster <- subset(coi_cluster,
                          cells = rownames(keep_clones))
    cl2_de <- run_de(coi_cluster, "condition", "Ri", "HD",
                     sprintf("clones_cluster_%s_Ri_v_HD", coi_cluster_id))
    
    coi_cluster_pb <- make_pseudobulk(coi_cluster, "condition")
    coi_cluster_de_pb <- FindAllMarkers(coi_cluster_pb, test.use = "DESeq2")
    write_csv(coi_cluster_de_pb,
              file.path(args$results,
                        sprintf("tables/clones_cluster_%s_Ri_v_HD_pseudobulk",
                                coi_cluster_id)))
    
    # Clone of interest, disease versus remission ----
    coi <- subset(seurat_obj, coi == "Ri01")
    
    coi_de <- run_de(coi, "Sample",
                    "Ri01_dis", "Ri01_5m", "Ri01_dis_v_5m")
    
    
}

# ----------------------------------------------------------------------------
main(args)
