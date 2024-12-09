library("tidyverse")
library("scRepertoire")

vdj_dir <- "data/raw/VDJ"
contig_annots <- list.files(vdj_dir,
                            recursive = TRUE,
                            full.names = TRUE,
                            pattern = "filtered_contig_annotations.csv") #"all_contig_annotations.csv")
names(contig_annots) <- basename(dirname(contig_annots))

samples <- read_delim("data/cellbender_input.txt", col_names = FALSE)[[1]]
contig_annots <- contig_annots[samples]

# Reformat to tab separated   
#    gsub("\\s+", ";",
#                readLines("data/cellbender_input.txt")) %>%
#    as_tibble() %>%
#    dplyr::filter(nchar(value) >= 1) %>%
#    tidyr::separate(value,
#                    into = c("sample", "cells", "droplets"),
#                    sep = ";",
#                    extra = "drop") %>%
#    write_delim("data/cellbender_input.txt", delim = "\t")
   
contig_list <- lapply(contig_annots, read_csv)

# Save the beta-chain cdr3
bchain_cdr3 <- lapply(names(contig_list), function(smp){
    contig_list[[smp]] %>%
        dplyr::filter(chain == "TRB") %>%
        dplyr::mutate(barcode = paste(smp, barcode, sep = "_")) %>%
        dplyr::select(barcode, cdr3)
})

# Save beta chain CDR3 ----

bchain_cdr3 <- dplyr::bind_rows(bchain_cdr3)
write_csv(bchain_cdr3,
          file = "data/processed/barcode_to_beta_chain_cdr3.csv")

# Combine clonotypes, save ----

combined_tcr <- combineTCR(contig_list, 
                           samples = names(contig_annots))
readr::write_rds(combined_tcr,
                 file = "data/processed/combined_tcr.rds")


# Plots of clone abundance ----
ri_samples <- grepl("^Ri", names(combined_tcr))
ri_nms <- names(combined_tcr)[ri_samples]
    
# Relative proportions of top Ri clones ----
#
#pdf("figures/scRepertoire/ri_top_10_clone_ppns.pdf", width = 10)
#p <- clonalCompare(combined_tcr, 
#              top.clones = 10, 
#              samples = ri_nms, 
#              cloneCall="aa", 
#              graph="alluvial") +
#    guides(fill=guide_legend(ncol = 1)) +
#    theme(axis.text = element_text(size = 14),
#          axis.title = element_text(size = 16),
#          legend.text=element_text(size = 10))
#print(p)
#dev.off()

# Scatter of proportions in Ri samples ----
#pdf(sprintf("figures/scRepertoire/clone_ppns_%s_%s.pdf", ri_nms[1], ri_nms[2]), width = 10)
#p <- clonalScatter(combined_tcr, 
#              cloneCall = "strict", 
#              x.axis = ri_nms[1],
#              y.axis = ri_nms[2],
#              dot.size = "total",
#              graph = "proportion") +
#    theme(axis.text = element_text(size = 14),
#          axis.title = element_text(size = 20),
#          axis.title.x = element_text(margin = margin(t = 10)),
#          axis.title.y = element_text(margin = margin(r = 10)))
#print(p)
#dev.off()

# Clonal homeostasis ----

pdf("figures/scRepertoire/clonal_homeostasis.pdf", width = 9)
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

pdf("figures/scRepertoire/clonal_rank_ppn.pdf", width = 9)
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

pdf("figures/scRepertoire/tcr_gene_usage.pdf", width = 8)
p <- 
vizGenes(combined_tcr, plot = "barplot") +
    theme(axis.text.x = element_text(size = 10))
print(p)
dev.off()    

