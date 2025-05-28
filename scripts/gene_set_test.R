library("edgeR")
library("limma")
library("tidyverse")


source(file.path(".", "scripts/markers_sp1.R"))
source(file.path(".", "scripts/funcs_custom_heatmaps.R"))

fig5_markers <- fig_5_markers()
all_markers <- unname(unlist(fig5_markers))
marker_df <- tibble(gene = unlist(fig5_markers),
                    category = rep(names(fig5_markers), lengths(fig5_markers)))

# make heatmap ----
make_heatmap <- function(dat,
                         col_group,
                         col_labs = NULL,
                         disp_min = -2.5,
                         disp_max = 2.5,
                         row_group = FALSE,
                         palette = viridis::viridis(100),
                         column_title_gp = gpar(fontsize = 10),
                         column_names_gp = gpar(fontsize = 10),
                         row_names_gp = gpar(fontsize = 3),
                         ...){
    if (! isTRUE(row_group) & ! isFALSE(row_group)) { 
        row_group <- cluster_within_group(t(dat), row_group)
    }
    
    
    Heatmap(dat,
            cluster_columns = cluster_within_group(dat, col_group),
            cluster_rows = row_group,
            show_column_dend = FALSE,
            show_row_dend = FALSE, 
            column_split = length(unique(col_group)),
            column_title_gp = column_title_gp,
            column_names_gp = column_names_gp, 
            row_names_gp = row_names_gp,
            column_labels = gsub("_.*", "", colnames(dat)),
            col = palette,
            heatmap_legend_param = list(title = "Scaled \nexpression"),
            ...)
    
}

# run edgeR ----
run_edger <- function(pb_fname, out_dir, min_total = 25,
                      min_reads = 2, min_samples = 3){
    
    pb_counts <- read_csv(pb_fname)
    rn <- pb_counts$gene
    pb_counts <- pb_counts %>% dplyr::select(-gene)
    cnd <- gsub(".*_(HD|Ri)$", "\\1", colnames(pb_counts))
    sample <- gsub("([^_]*)_.*", "\\1", colnames(pb_counts))
    pb <- as.matrix(pb_counts, dimnames = list(rn, colnames(pb_counts)))
    rownames(pb) <- rn
    
    samples <- data.frame(sample=sample,
                         condition=cnd)
    design <- model.matrix(~0 + condition, data = samples)
    y <- DGEList(counts = pb, samples = samples)
 
    keep_genes <- rowSums(pb) >= min_total &
        rowSums(pb >= min_reads) >= min_samples
    
    y <- y[keep_genes, , keep=FALSE]
    print(table(all_markers %in% rownames(y)))
    
    y <- normLibSizes(y)
    y <- estimateDisp(y, design, robust=TRUE)
    plotBCV(y)
    
    fit <- glmQLFit(y, design, robust=TRUE)
    plotQLDisp(fit)
    
    qlf <- glmQLFTest(fit)
    de_res <- topTags(qlf, n = nrow(pb))  
    sig <- de_res$table %>%
        dplyr::filter(FDR <= 0.05)
                         
    # MDS
    cluster <- as.factor(y$samples$sample)
    plotMDS(y, pch=16, col=c(2:10)[cluster], main="MDS")
    legend("topleft", legend=paste0(levels(cluster)),
           pch=16, col=2:8, cex=0.8)
   
     logcpm <- cpm(y[intersect(rownames(sig), all_markers), ], log=TRUE)
     log_cpm_cats <- marker_df[match(rownames(logcpm), marker_df$gene), ]$category
     
     write_csv(as_tibble(logcpm, rownames = "gene"),
               file = file.path(out_dir, "edgeR_log_cpm.csv"))
     write_csv(as_tibble(de_res$table, rownames = "gene"),
                         file = file.path(out_dir, "edgeR_de_res.csv"))
     
     h <- make_heatmap(logcpm,
                       col_group = cnd,
                       row_names_gp = gpar(fontsize = 6),
                       row_group = log_cpm_cats,
                       row_split = length(unique(log_cpm_cats)),
                       row_title_gp = gpar(fontsize = 7.4))
     
 
     pdf(file.path(out_dir, "edgeR_heatmap_cpm_sig_genes_by_cat.pdf"),
         height = 12)
     print(h)
     dev.off()
     
     # Limma
     dge <- DGEList(counts = pb[rownames(y),], samples = samples)
     dge <- calcNormFactors(dge, method = "TMMwsp")
     
     vfit2 <- voomLmFit(dge, design=design)#, sample.weights=TRUE)
     vfit2 <- eBayes(vfit2)
     tt <- topTable(vfit2, coef=2, sort.by="P", n = nrow(vfit2))
     
        
}

# Analysis ----



base_dir <- "reactive_cl1/"
cd3 <- file.path(base_dir, "reactive_Ri_v_HD_pseudo_cdr3")
clone <- file.path(base_dir, "reactive_Ri_v_HD_pseudo_clone")
smp <- file.path(base_dir, "reactive_Ri_v_HD_pseudo_sample")

run_edger(sprintf("%s/summed_counts.csv", cd3), cd3)
run_edger(sprintf("%s/summed_counts.csv", clone), clone)
run_edger(sprintf("%s/summed_counts.csv", smp), smp)

     
