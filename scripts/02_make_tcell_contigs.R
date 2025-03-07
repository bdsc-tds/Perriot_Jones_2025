library("tidyverse")
library("scRepertoire")

vdj_dir <- "data/raw/VDJ"
figs <- "figures/scRepertoire/"
data_out <- "data/processed/"
    
contig_annots <- list.files(vdj_dir,
                            recursive = TRUE,
                            full.names = TRUE,
                            pattern = "all_contig_annotations.csv")
names(contig_annots) <- basename(dirname(contig_annots))

samples <- read_delim("data/cellbender_input.txt", col_names = FALSE)[[1]]
all(samples %in% names(contig_annots))
contig_annots <- contig_annots[samples]
contig_list <- lapply(contig_annots, read_csv)

# Save the beta-chain CDR3 ----
bchain_cdr3 <- lapply(names(contig_list), function(smp){
    contig_list[[smp]] %>%
        dplyr::filter(chain == "TRB") %>%
        dplyr::mutate(barcode = paste(smp, barcode, sep = "_")) %>%
        dplyr::select(barcode, cdr3)
})

bchain_cdr3 <- dplyr::bind_rows(bchain_cdr3)
write_csv(bchain_cdr3,
          file = file.path(data_out, "barcode_to_beta_chain_cdr3.csv"))

# Combine clonotypes, save ----

combined_tcr <- combineTCR(contig_list, 
                           samples = names(contig_annots))
readr::write_rds(combined_tcr,
                 file = file.path(data_out, "combined_tcr.rds"))


# Plots of clone abundance ----

# Clonal homeostasis ----

pdf(file.path(figs, "clonal_homeostasis.pdf"), width = 9)
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

pdf(file.path(figs, "clonal_rank_ppn.pdf"), width = 9)
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

pdf(file.path(figs, "tcr_gene_usage.pdf"), width = 8)
p <- vizGenes(combined_tcr, plot = "barplot") +
    theme(axis.text.x = element_text(size = 10))
print(p)
dev.off()    

