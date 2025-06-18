# filter_tcr_genes ----
filter_tcr_genes <- function(markers){
    markers %>%
        dplyr::filter(! grepl("^TR[ABG][VDJ]", gene))
} 


# Pairwise differential expression ----
pairwise_de <- function(seurat_obj,
                        results,
                        pairwise = de_expanded(),
                        test_use = "DESeq2",
                        p_cutoff = 0.05,
                        n_top = 50,
                        min_pct = 0.1, 
                        remove_tcr = TRUE){
    
    # Create filenames using options ----
    if (isTRUE(remove_tcr)){
        filter_st <- "_no_tcrs_"
    } else {
        filter_st <- ""
    }
    heatmap_template <- file.path(results, "heatmap_pairwise_%stop_%s%s.pdf")
     
    if (test_use == "DESeq2"){ 
        res_template <- file.path(results, "de_%s_pseudobulk.csv")
        heatmap_fname <- sprintf(heatmap_template, "pseudobulk",
                                 n_top, filter_st)
    } else {
        res_template <- file.path(results, "de_one_v_others_%s.csv")
        heatmap_fname <- sprintf(heatmap_template, "one_v_others",
                                 n_top, filter_st)
    }
    
    # Run differential expression ----    
    markers <- lapply(seq_len(nrow(pairwise)), function(i){
        markers <- FindMarkers(seurat_obj,
                               ident.1 = pairwise[[i, "ident_1"]],
                               ident.2 = pairwise[[i, "ident_2"]], 
                               test.use = test_use) %>%
            tibble::as_tibble(rownames = "gene") %>%
            dplyr::mutate(condition = pairwise[[i, "name"]]) %>%
            dplyr::arrange(desc(abs(avg_log2FC))) %>%
            dplyr::relocate(condition, gene)
        
        write_csv(markers, sprintf(res_template, pairwise[[i, "name"]]))
        markers %>%
            dplyr::filter(p_val_adj <= p_cutoff)
    })
    markers <- dplyr::bind_rows(markers)
    
    if (isTRUE(remove_tcr)){ markers <- filter_tcr_genes(markers) } 
    
    pairwise_de_heatmap(seurat_obj, markers, heatmap_fname, n_top)
    invisible(markers)
}

# Heatmap of pairwise differential expression ----
pairwise_de_heatmap <- function(seurat_obj, markers, outfname, n_top){
    heatmap_markers <- markers %>%
        dplyr::group_by(condition) %>%
        dplyr::slice_head(n = n_top)
    
    pdf(outfname)
    p <- DoHeatmap(seurat_obj,
                   size = 4,
                   angle = 90,
                   features = unique(heatmap_markers$gene)) +
        theme(axis.text = element_text(size = 2))
    print(p)
    dev.off()
}