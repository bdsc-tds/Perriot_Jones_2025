# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ComplexHeatmap")
library("DESeq2")
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
source(file.path(args$workdir, "scripts/funcs_custom_heatmaps.R"))
source(file.path(args$workdir, "scripts/clones_of_interest.R"))
# ----------------------------------------------------------------------------

# make_pseudobulk - wrapper with corrections for when a category has no counts
# for any gene ----
make_pseudobulk <- function(obj, idents, group_by = c("Sample", "beta_aa")){
    gp_by = unique(c(group_by, idents))
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
make_volcano <- function(de, res_dir){
    pdf(file.path(res_dir, "volcano.pdf"))
    p <- volc_mod(de,
                  xlab = expression("Average " * log[2] * "FC"),
                  ylab = expression(-log[10]*"(adjusted p-value)"),
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

    de <- de %>% filter(! is.na(p_val_adj))
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
     
    # Pseudobulk heatmap ----
    if (isTRUE(pseudobulk)) { obj_subs <- make_pseudobulk(obj_subs, idents) } 
    
    pdf(gsub(".pdf", "_pseudobulk_clustered.pdf", out_fname),
        width = wd, height = ht)
    h <- heatmap_w_labs(obj_subs,
                        col_labs = structure(levels(Idents(obj_subs)),
                                             names = levels(Idents(obj_subs))), 
                        col_group = idents, 
                        row_group = FALSE)
    print(h)
    dev.off()
}

# Run differential expression analyses ----
run_diff_expr <- function(obj, results_dir, idents, ident_1, ident_2,
                          name, max_n = 100,
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
    make_volcano(de, res_dir)
        
    # CSV of gene expression
    write_top_genes(de, obj, res_dir, pval_min, fc_cutoff = 0.1)

    tryCatch({dev.off()}, error = function(e){return()})
    
    # Heatmap
    heatmap_out <- file.path(res_dir, "heatmap.pdf")

    make_heatmaps(obj = obj, de = de, idents = idents, out_fname = heatmap_out,
                  max_n = max_n, pval_min = pval_min, wd = wd, ht = ht,
                  group.by = idents)
    
    return(de)
}

# run_diff_expr_pb ----
run_diff_expr_pb <- function(obj, results_dir, idents, ident_1, ident_2,
                             group_by, name, max_n = 100, test_use = "DESeq2",
                             pval_min = 0.05, wd = 7, ht = 7, ...){
    res_dir <- file.path(results_dir, name)
    if (! file.exists(res_dir)) { dir.create(res_dir, recursive = TRUE) }
    
    pseudo <- make_pseudobulk(obj, idents, group_by)
    
    # Write summed and average counts ----
    pseudo_counts <- as.matrix(LayerData(pseudo, "counts")) %>%
        as_tibble(rownames = "gene")
    write_csv(pseudo_counts, file.path(res_dir, "summed_counts.csv"))

    # Write DESeq2 equivalent of cpm ----
    
    
    
    
    
    av_expr <- as.matrix(AverageExpression(obj,
                               group.by = unique(c(idents, group_by)),
                               layer = "counts")$RNA) %>%
        as_tibble(rownames = "gene")
    write_csv(av_expr, file.path(res_dir, "average_counts.csv"))
    
    # Run diff expr ----
    Idents(pseudo) <- pseudo[[idents]][, 1]
    de <- FindMarkers(pseudo, ident.1 = ident_1, ident.2 = ident_2, 
                      test.use = test_use, ...) %>%
        dplyr::as_tibble(rownames = "gene") %>%
        dplyr::mutate("comparison" = name) %>%
        relocate(gene, comparison)
    
    # CSV of results
    write_csv(de, file.path(res_dir, "diff_expr_pseudo_full.csv"))
    
    de <- filter_tcr_genes(de)
    write_csv(de, file.path(res_dir, "diff_expr_pseudo_rm_tcr.csv"))
    
    make_heatmaps(pseudo, de, idents,
                  out_fname = file.path(res_dir, "heatmap_pseudobulk_de.pdf"),
                  max_n = max_n,
                  pseudobulk = FALSE)
    
    return(de)
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

fig_5f_all <- function(){
    markers <- c("HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", 
                 "IL2RB", "CD69", "CD74", "TFRC", "CD44", "TNFRSF9", "VCAM1", 
                 "ICAM1", "TNF", "IFNG", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK", 
                 "FASLG", "LAMP1", "NKG7", "CCL4", "CCL3", "CSF1", "IFIT1", "RSAD2", 
                 "IFIT3", "MX1", "OAS3", "IFI6", "OAS1", "ISG15", "IFI44L", "IFIT2", 
                 "IKZF2", "KLRC1", "KLRC2", "KLRC3", "KLRK1", "KIR2DL1", "KIR2DL2", 
                 "KIR2DL3", "KIR2DL4", "KIR3DL1", "KIR3DL2", "KIR3DL3", "IL7R", 
                 "ZAP70", "LCK", "NFATC1", "NFKB1", "NFKB2", "FOS", "JUN", "RELB", 
                 "EGR1", "EGR2", "EGR3", "JUND", "FOSB", "PTPN6", "PTPN11", "INPP5D", 
                 "CBL", "SHC1", "DOK2", "RASA1", "ITCH", "CBLB", "RC3H1", "RC3H2", 
                 "DGKE", "DGKA", "ZFP36L1", "NFKBIA", "SH2B3", "EMPP1", "LRRK1", 
                 "IRF8", "PLCG2", "ABCB9", "PRKCE", "DAPK2", "ID2", "ADAM9", "EPAS1"
    )
    return(markers)
}


make_violins <- function(seurat_obj, args, out_dir = ".", markers = NULL){
    if (is.null(markers)) {
        source(file.path(args$workdir, "scripts/markers_sp1.R"))
        markers <- fig_5_markers()
        
    }
    
    # Violin plots using log counts
    x <- lapply(names(markers), function(nm) {
        pdf(file.path(out_dir, sprintf("%s_violin.pdf", gsub(" ", "_", nm))),
            width = 12, height = 20)
        h <- VlnPlot(seurat_obj, features = markers[[nm]],
                     group.by = "Sample", ncol = 2)
        print(h)
        dev.off() 
    })
    
    
    # Violin plots using raw counts
    x <- lapply(names(markers), function(nm) {
        pdf(file.path(out_dir,
                      sprintf("%s_violin_counts.pdf", gsub(" ", "_", nm))),
            width = 12, height = 20)
        h <- VlnPlot(seurat_obj, features = markers[[nm]],
                     group.by = "Sample", ncol = 2, layer = "counts")
        print(h)
        dev.off() 
    })
    
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
    
    # Write table of counts per cluster ---- 
    cluster_counts_table(seurat_obj[[]] %>%
                             dplyr::filter(seurat_clusters == coi_cluster_id),
                      file.path(args$results,
                                sprintf("tables/clone_counts_cl_%s.csv",
                                        coi_cluster_id)))
    
    run_de <- purrr::partial(run_diff_expr, results = args$results, ...)
    run_pseudo <- purrr::partial(run_diff_expr_pb, results = args$results, ...)
    
    # Subset to clones where reactivity information is known ----
    clones <- subset(seurat_obj, reactive %in% c(TRUE, FALSE))
    write_rds(clones, file.path(dirname(args$seurat), "reactive_clones.rds"))
    
    # Write table of reactive counts ----
    clones_counts_table(clones[[]],
                        file.path(args$results,
                                  "tables/reactive_clone_counts.csv"))
    
    # Write a table of cluster counts for clones from cluster of interest ----
    cl_clones <- seurat_obj[[]] %>%
        dplyr::filter(seurat_clusters == coi_cluster_id) %>%
        dplyr::select(Sample, beta_aa) %>%
        unique()
    cl_clusters <- seurat_obj[[]] %>%
        dplyr::inner_join(cl_clones, relationship = "many-to-one") %>%
        dplyr::group_by(Sample, beta_aa, tcr_name) %>%
        dplyr::mutate(n_clone = n()) %>%
        dplyr::group_by(Sample, beta_aa, tcr_name, seurat_clusters, n_clone) %>%
        dplyr::summarise(n_clone_cluster = n()) %>%
        dplyr::arrange(desc(n_clone), Sample, beta_aa)
    
    write_csv(cl_clusters,
              file.path(args$results,
                        sprintf("tables/clusters_with_cells_in_cl_%s.csv",
                                coi_cluster_id)))
    
    write_csv(cl_clusters %>%
                   tidyr::pivot_wider(names_from = seurat_clusters,
                                      values_from = n_clone_cluster),
              file.path(args$results,
                        sprintf("tables/clusters_with_cells_in_cl_%s_wide.csv",
                                coi_cluster_id)))
    
    
    # Remove followup sample, adjust condition ---- 
    clones_wo_5m <- subset(clones, Sample != "Ri01_5m")
    
    clones_wo_5m[[]] <- clones_wo_5m[[]] %>%
        dplyr::mutate(condition = case_when(condition == "Ri01_dis" ~ "Ri",
                                            TRUE ~ condition),
                      reactive_lab = case_when(reactive == "TRUE" ~ "Reactive",
                                               reactive == "FALSE" ~ "Non-reactive"))
    
    clones_wo_5m[[]] <- clones_wo_5m[[]] %>%
        dplyr::mutate(rx_by_cnd = paste(condition, reactive, sep = "_"))
    Idents(clones_wo_5m) <- "rx_by_cnd"
    
    #_______________________________________________________
    # Violin plots and heatmaps for cluster 1 reactive 
    
    clones_cl1 <- subset(clones_wo_5m, seurat_clusters == "1" & reactive == "TRUE")
    
    markers <- fig_5f_all()
    
    rx_cl1_markers <- subset(clones_cl1, features = markers)
    rx_cl1_markers <- ScaleData(rx_cl1_markers,
                                features = Features(rx_cl1_markers))
    
    pdf("rx_cl1.pdf")
    h <- DoHeatmap(rx_cl1_markers, features = Features(rx_cl1_markers))
    print(h)
    dev.off()
    
    
    pdf("rx_cl1_by_donor.pdf", width = 12)
    h <- heatmap_w_labs(rx_cl1_markers, col_group = "Sample")
    print(h)
    dev.off()
    
    # Make violins ----
    
    # For just cluster 1
    source(file.path(args$workdir, "scripts/markers_sp1.R"))
    fig5_markers <- fig_5_markers()
    
    make_violins(rx_cl1_markers, args, outdir = ".", markers = fig5_markers)
    
    # For all reactive 
    all_rx <- subset(clones_wo_5m,
                     reactive == "TRUE",
                     features = unlist(markers))
    all_rx <- ScaleData(all_rx,
                        features = Features(all_rx))
    
    make_violins(all_rx, args,
                 out_dir = "results/violins/all_reactive/",
                 markers = fig5_markers)
    #_______________________________________________________
    # Compare cluster 1 versus others 
    all_rx[[]] <- all_rx[[]] %>%
        dplyr::mutate(is_coi = ifelse(seurat_clusters == 1, TRUE, FALSE))
    Idents(all_rx) <- all_rx$is_coi
    cell_de <- FindAllMarkers(all_rx) %>%
        as_tibble() %>%
        dplyr::rename(cluster_1 = cluster)
    write_csv(cell_de,
              file.path("results/reactive_cl1_v_others/cell_level_diff_expr.csv"))
    
    all_rx_pb <- AggregateExpression(all_rx,
                                     group.by = c("Sample", "is_coi"),
                                     return.seurat = TRUE)
    Idents(all_rx_pb) <- all_rx_pb$is_coi
    sample_pb <-  FindAllMarkers(all_rx_pb, test.use = "DESeq2") # No DE genes
    write_csv(sample_pb,
              file.path("results/reactive_cl1_v_others/sample_pseudobulk_diff_expr.csv"))
    
    
    all_rx_cdr3_pb <- AggregateExpression(all_rx,
                                          group.by = c("TCR2", "is_coi"),
                                          return.seurat = TRUE)
    
    Idents(all_rx_cdr3_pb) <- all_rx_cdr3_pb$is_coi
    cdr3_pb <-  FindAllMarkers(all_rx_pb, test.use = "DESeq2") # No DE genes
    write_csv(cdr3_pb,
              file.path("results/reactive_cl1_v_others/cdr3_pseudobulk_diff_expr.csv"))
    
    # Barplot reactive all clusters
    pdf("results/reactive_cl1_v_others/barplot_all_clusters.pdf")
    all_rx[[]] %>%
        ggplot(aes(x = seurat_clusters, fill = Sample)) +
        geom_bar(position = position_fill()) + theme_bw()
    dev.off()
    
    # Barplot reactive cl1 v others
    pdf("results/reactive_cl1_v_others/barplot_cl1_v_others.pdf")
    all_rx[[]] %>%
        ggplot(aes(x = is_coi, fill = Sample)) +
        geom_bar(position = position_fill()) + theme_bw()
    dev.off()
    
    pdf("results/reactive_cl1_v_others/barplot_cl1_v_others_by_condition.pdf")
    all_rx[[]] %>%
        mutate(cl_cond = paste(condition, is_coi, sep = "_")) %>%
        ggplot(aes(x = cl_cond, fill = Sample)) +
        geom_bar(position = position_fill()) + theme_bw() +
        labs(x = "Cluster 1 versus other clusters") 
    dev.off()
    
    #_______________________________________________________
    
    # Subset to just reactive clones, test HD v Ri ----
    rx <- subset(clones_cl1, reactive == "TRUE")
    cnd_de <- run_de(rx, "condition", "Ri", "HD", "reactive_Ri_v_HD")
    
    # Pseudobulk at the sample level
    cnd_de_pb_sample <- run_pseudo(rx, "condition", "Ri", "HD", 
                                   group_by = "Sample",
                                   name = "reactive_Ri_v_HD_pseudo_sample")
    
    # Pseudobulk by complementarity determining region CDR3
    cnd_de_pb_cdr3 <- run_pseudo(rx, "condition", "Ri", "HD", 
                                 group_by = c("Sample", "CTaa2"),
                                 name = "reactive_Ri_v_HD_pseudo_cdr3")
    
    # Pseudobulk by CDR3 and beta chains
    cnd_de_pb_beta_aa <- run_pseudo(rx, "condition", "Ri", "HD", 
                                    group_by = c("Sample", "beta_aa"),
                                    name = "reactive_Ri_v_HD_pseudo_clone")
    
    
    # Pseudobulk by complementarity determining region CDR3 with Wilcoxon
    cnd_de_pb_cdr3 <- run_pseudo(rx, "condition", "Ri", "HD", 
                                 group_by = c("Sample", "CTaa2"),
                                 test_use = "wilcox",
                                 name = "reactive_Ri_v_HD_pseudo_cdr3_wilcox")
    
    # Pseudobulk by clone with Wilcoxon
    cnd_de_pb_clone <- run_pseudo(rx, "condition", "Ri", "HD", 
                                  group_by = c("Sample", "beta_aa"),
                                  test_use = "wilcox",
                                  name = "reactive_Ri_v_HD_pseudo_clone_wilcox")
    
    
    # --------------------------
    # Delme?
    

    # One v all pseudobulk - aggregate samples + beta_aa ---- 
    
    # This includes filtering
    pb_beta_aa <- make_pseudobulk(clones_wo_5m, "rx_by_cnd")
    
    pb_beta_aa_de <- FindAllMarkers(pb_beta_aa,
                                    test.use = "DESeq2") %>%
        as_tibble(rownames = "gene")
    
    write_csv(pb_beta_aa_de,
              file.path(args$results,
                        "condition_by_reactivity/de_full_pseudobulk_beta_aa.csv"))
    
    
    # Pseudobulk samples ----
    pb_sample <- make_pseudobulk(clones_wo_5m, "rx_by_cnd", group_by = "Sample")
    
    pb_sample_de <- FindAllMarkers(pb_sample,
                                   test.use = "DESeq2") %>%
        as_tibble(rownames = "gene")
    
    write_csv(pb_sample_de,
              file.path(args$results,
                        "condition_by_reactivity/de_full_pseudobulk_sample.csv"))
    
    # Pseudobulk CDR3 ----
    pb_cdr3 <- make_pseudobulk(clones_wo_5m, "rx_by_cnd",
                               group_by = c("Sample", "CTaa2"))
    
    pb_cdr3_de <- FindAllMarkers(pb_cdr3,
                                 test.use = "DESeq2") %>%
        as_tibble(rownames = "gene")
    
    write_csv(pb_sample_de,
              file.path(args$results,
                        "condition_by_reactivity/de_full_pseudobulk_cdr3.csv"))
}

# ----------------------------------------------------------------------------
main(args)

# Add TCR annotation to heatmap  
# 
# so <- seurat_obj[[]] %>%
#     dplyr::select(beta_aa, matches("trb"), cdr3_aa_beta) %>%
#     unique() %>%
#     dplyr::mutate(beta_aa = gsub("_", "-", beta_aa))
# 
# clone_pseudo[[]] <- clone_pseudo[[]] %>%
#     dplyr::left_join(so) 
# 
# 
# column_ha = HeatmapAnnotation(TRBV = clone_pseudo$trbv,
#                               TRBJ = clone_pseudo$trbj)
# 
# h <- heatmap_w_labs(clone_pseudo,
#                     col_labs = structure(levels(Idents(clone_pseudo)),
#                                          names = levels(Idents(clone_pseudo))),
#                     col_group = idents,
#                     row_group = FALSE, show_column_names = TRUE)
# 
# 
# column_ha = HeatmapAnnotation(CDR3 = anno_text(clone_pseudo$cdr3_aa_beta,
#                                                rot = 90,
#                                                location = 0,
#                                                just = "left",
#                                                gp = gpar(fontsize = 8))