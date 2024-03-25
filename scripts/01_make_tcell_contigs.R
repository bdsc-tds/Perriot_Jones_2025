library("tidyverse")
library("scRepertoire")

vdj_dir <- "data/raw/VDJ"
contig_annots <- list.files(vdj_dir, recursive = TRUE, full.names = TRUE)
names(contig_annots) <- basename(dirname(contig_annots))

contig_list <- lapply(contig_annots, read_csv)

# Save the beta-chain cdr3
bchain_cdr3 <- lapply(names(contig_list), function(smp){
    contig_list[[smp]] %>%
        dplyr::filter(chain == "TRB") %>%
        dplyr::mutate(barcode = paste(smp, barcode, sep = "_")) %>%
        dplyr::select(barcode, cdr3)
})

bchain_cdr3 <- dplyr::bind_rows(bchain_cdr3)
write_csv(bchain_cdr3, file = "results/barcode_to_beta_chain_cdr3.csv")

combined_tcr <- combineTCR(contig_list, 
                           samples = names(contig_annots))


# Plots of clone abundance ----
ri_samples <- grepl("^Ri", names(combined_tcr))
ri_nms <- names(combined_tcr)[ri_samples]
    
# Relative proportions of top Ri clones ----

pdf("figures/ri_top_10_clone_ppns.pdf", width = 10)
p <- clonalCompare(combined_tcr, 
              top.clones = 10, 
              samples = ri_nms, 
              cloneCall="aa", 
              graph="alluvial") +
    guides(fill=guide_legend(ncol = 1)) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text=element_text(size = 10))
print(p)
dev.off()


# Scatter of proportions in Ri samples ----
pdf(sprintf("figures/clone_ppns_%s_%s.pdf", ri_nms[1], ri_nms[2]), width = 10)
p <- clonalScatter(combined_tcr, 
              cloneCall = "strict", 
              x.axis = ri_nms[1],
              y.axis = ri_nms[2],
              dot.size = "total",
              graph = "proportion") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 20),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)))
print(p)
dev.off()

# Clonal homeostasis ----

pdf("figures/clonal_homeostasis.pdf", width = 9)
p <- 
clonalHomeostasis(combined_tcr) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10))) +
    labs(x = NULL)
print(p)
dev.off()

pdf("figures/clonal_rank_ppn.pdf", width = 9)
p <- 
clonalProportion(combined_tcr, cloneCall = "strict") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 15)),
          axis.title.y = element_text(margin = margin(r = 15)))
print(p)
dev.off()


# TCR gene usage ----

pdf("figures/tcr_gene_usage.pdf", width = 8)
p <- 
vizGenes(combined_tcr, plot = "barplot") +
    theme(axis.text.x = element_text(size = 10))
print(p)
dev.off()    

# Write out combined clonotypes ----
readr::write_rds(combined_tcr, file = "data/processed/combined_tcr.rds")

# Experimental ----
#clonalQuant(combined_tcr,#[ri_samples], 
#            cloneCall="strict", 
#            scale = TRUE)

#clonalAbundance(combined_tcr[ri_samples], 
#                cloneCall = "strict", 
#                scale = FALSE)

#clonalLength(combined_tcr[ri_samples], 
#             cloneCall="aa", 
#             chain = "both") 

#vizGenes(combined_tcr) +
#    theme(axis.text.y = element_text(size = 16)) 