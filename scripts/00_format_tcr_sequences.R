# Libraries and setup ----

library("janitor")
library("readxl")
library("tidyverse")

alpha_updates <- c("TRAV29"="TRAV29/DV5")

# only in combination in the md 
#TRAV12-2.TRAJ23.CAVNQGGKLIF_TRBV6-2.TRBJ2-2.CASSFEGAGELFF;TRAV21.TRAJ31.CAVRVNNARLMF_TRBV6-2.TRBJ2-2.CASSFEGAGELFF

raw_tcr <- "data/raw/20250114_TCR_list_Helen.xlsx"

t <- readxl::read_excel(raw_tcr, sheet = 1)  %>%
    janitor::remove_empty(which = c("rows", "cols")) %>%
    janitor::clean_names() %>%
    dplyr::mutate(sample = gsub("_.*", "", name)) %>%
    dplyr::relocate(sample) %>%
    
    # Rename the CDR3 columns 
    dplyr::rename("cdr3_aa_alpha" = cdr3_aaseq_3,
                  "cdr3_aa_beta" = cdr3_aaseq_6) %>%
    
    dplyr::mutate(alpha_vj_aa = paste(trav, traj, cdr3_aa_alpha, sep = "."),
                  beta_vj_aa = paste(trbv, trbj, cdr3_aa_beta, sep = "."),
                  vj_aa = paste(alpha_vj_aa, beta_vj_aa, sep = "_")) 

# Check all found ----

md <- read_csv("data/processed/cd8_and_tcr/integrated_seurat.csv.gz") %>%
    dplyr::select(Sample, matches("CT|TCR"), vj_aa) %>%
    dplyr::mutate(beta_vj_aa = gsub(".*_", "", vj_aa)) 

missing <- t %>% filter(! t$vj_aa %in% md$vj_aa)
semi_join(md, missing, by = "beta_vj_aa") %>% unique()

# Patch differences -----

t <- t %>%
    dplyr::rename(tcr_name = name) %>%
    dplyr::mutate(across(c(trav, alpha_vj_aa, vj_aa),
                         ~str_replace_all(.x, alpha_updates)))

# Write processed table ----
write_csv(t, file = "data/processed/tcr_seqs.csv")

