# Libraries and setup ----

library("janitor")
library("readxl")
library("tidyverse")

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
    dplyr::select(Cell, Sample, matches("CT|TCR"), vj_aa) %>%
    dplyr::mutate(beta_vj_aa = gsub(".*_", "", vj_aa))

missing <- t %>% filter(! t$vj_aa %in% md$vj_aa)
semi_join(md, missing, by = "beta_vj_aa")

# Write processed table ----
m