# filter_tcr_genes ----
filter_tcr_genes <- function(markers){
    markers %>%
        dplyr::filter(! grepl("^TR[ABG][VDJ]", gene))
} 


# make_pseudobulk - wrapper with corrections for when a category has no counts
# for any gene ----
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


# run_diff_expr_pb ----
run_diff_expr_pb <- function(obj, results_dir, idents, ident_1, ident_2,
                             name, markers = functional_markers(), max_n = 100,
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
