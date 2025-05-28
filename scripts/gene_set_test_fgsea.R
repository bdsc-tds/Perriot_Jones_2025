library("fgsea")
library("tidyverse")

source(file.path(".", "scripts/markers_sp1.R"))

fig5_markers <- fig_5_markers()
all_markers <- unname(unlist(fig5_markers))
marker_df <- tibble(gene = unlist(fig5_markers),
                    category = rep(names(fig5_markers), lengths(fig5_markers)))


sample_de <- "reactive_cl1/reactive_Ri_v_HD_pseudo_sample/diff_expr_pseudo_rm_tcr.csv"
de_res <- read_csv(sample_de) %>%
    dplyr::arrange(p_val)

de_ranks <- structure(de_res$p_val, names = de_res$gene)
de_ranks <- de_ranks[! is.na(de_ranks)]

# Need a p-value cutoff for this
#fc_ranks <- structure(de_res$avg_log2FC, names = de_res$gene)
#fc_ranks <- fc_ranks[! is.na(fc_ranks)]


fgseaRes <- fgsea(pathways = fig5_markers, 
                  stats    = de_ranks,
                  scoreType = "pos") # using pos as giving P-values not logFC


plotGseaTable(fig5_markers, de_ranks, fgseaRes, gseaParam=0.5)
