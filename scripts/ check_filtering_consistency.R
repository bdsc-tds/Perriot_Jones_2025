library("Seurat")
library("tidyverse")

# Functions -----

# make_wide_md ----
make_long_md <- function(md){
    md %>%
        dplyr::select(Sample, seurat_clusters, Cell, TCR2, CTaa2, beta_aa) %>%
        # TCR beta must be present
        dplyr::filter(! is.na(TCR2)) %>%
        dplyr::mutate(across(c(TCR2, CTaa2), ~strsplit(.x, ";"))) %>%
        tidyr::unnest(c(TCR2, CTaa2)) %>%
        tidyr::separate_wider_delim(cols = c(TCR2),
                                    delim = '.', names_sep = "_") %>%
        dplyr::rename(trbv = TCR2_1, 
                      trbj = TCR2_3,
                      cdr3_aa_beta = CTaa2) 
}


# ---------

workdir='/data/PRTNR/CHUV/DIR/rgottar1/encephalitis_tcr/Perriot_Jones/'
if (! file.exists(workdir)) { workdir <- "." }

selected_clones <- file.path(workdir, "data/processed/merged_clone_sets.csv")
clones <- read_csv(selected_clones)

# CD8 or tcr filtered ----
cd8_or_tcr <- file.path(workdir,
                        "data/processed/cd8_or_tcr/integrated_seurat.csv.gz")
cd8_or_tcr <- read_csv(cd8_or_tcr) 

cd8_or_tcr_long <- make_long_md(cd8_or_tcr) 

x <- cd8_or_tcr_long %>% dplyr::inner_join(clones)

# Check if single cell has two different reactive clones
x %>%
    dplyr::group_by(Cell) %>%
    dplyr::filter(n_distinct(tcr_name) > 1) # 13 are shared

cd8_or_tcr_clusters <- x %>%
    dplyr::group_by(tcr_name) %>%
    dplyr::mutate(n_clone = n()) %>%
    dplyr::group_by(tcr_name, seurat_clusters, n_clone) %>%
    dplyr::summarise(n_cluster = n()) %>%
    dplyr::mutate(pct_clone_cluster = n_cluster/n_clone*100) %>%
    dplyr::arrange(tcr_name, desc(pct_clone_cluster))

cd8_or_tcr_top <- cd8_or_tcr_clusters %>%
    dplyr::group_by(tcr_name) %>%
    dplyr::slice_head(n = 1)


# CD8 and tcr filtered ----
cd8_and_tcr <- paste0("/data/PRTNR/CHUV/DIR/rgottar1/encephalitis_tcr/",
                      "Perriot_Jones_v1/archive/Feb_2025/data/",
                      "processed/cd8_and_tcr/integrated_seurat.csv.gz")
cd8_and_tcr <- read_csv(cd8_and_tcr) 

cd8_and_tcr_long <- make_long_md(cd8_and_tcr) 

x <- cd8_and_tcr_long %>%
    dplyr::inner_join(clones)

# Check if single cell has two different reactive clones
x %>%
    dplyr::group_by(Cell) %>%
    dplyr::filter(n_distinct(tcr_name) > 1) # Yes 8, cells

cd8_and_tcr_clusters <- x %>%
    dplyr::group_by(tcr_name) %>%
    dplyr::mutate(n_clone = n()) %>%
    dplyr::group_by(tcr_name, seurat_clusters, n_clone) %>%
    dplyr::summarise(n_cluster = n()) %>%
    dplyr::mutate(pct_clone_cluster = n_cluster/n_clone*100) %>%
    dplyr::arrange(tcr_name, desc(pct_clone_cluster))
    
cd8_and_tcr_top <- cd8_and_tcr_clusters %>%
    dplyr::group_by(tcr_name) %>%
    dplyr::slice_head(n = 1)
    
# ----
one_tcr_beta <- file.path(workdir, "/data/processed/one_tcr_beta/integrated_seurat.csv.gz")
one_tcr_beta <- read_csv(one_tcr_beta)

bb_long <- make_long_md(one_tcr_beta) 
bb <- bb_long %>%
    dplyr::inner_join(clones)

# Check that all clones are in bb
clones %>% dplyr::anti_join(bb)

# Check if single cell has two different reactive clones
bb %>%
    dplyr::group_by(Cell) %>%
    dplyr::filter(n_distinct(tcr_name) > 1) # No, by definition

bb_clusters <- bb %>%
    dplyr::group_by(tcr_name) %>%
    dplyr::mutate(n_clone = n()) %>%
    dplyr::group_by(tcr_name, seurat_clusters, n_clone) %>%
    dplyr::summarise(n_cluster = n()) %>%
    dplyr::mutate(pct_clone_cluster = n_cluster/n_clone*100) %>%
    dplyr::arrange(tcr_name, desc(pct_clone_cluster))

bb_top <- bb_clusters %>%
    dplyr::group_by(tcr_name) %>%
    dplyr::slice_head(n = 1)

# ---------

colnames(bb_top) <- paste0(colnames(bb_top), "_one_beta")
bb_top <- bb_top %>% dplyr::rename(tcr_name = tcr_name_one_beta)

x <- cd8_or_tcr_top %>%
    
    dplyr::left_join(cd8_and_tcr_top,
                     by = "tcr_name",
                     suffix = c("_cd8_or_tcrb", "_cd8_and_tcr")) %>%
    
    dplyr::left_join(bb_top,
                     by = "tcr_name") %>%
 
    dplyr::relocate(tcr_name,
                    matches("n_clone"),
                    matches("n_cluster"),
                    matches("pct_clone_cluster")) 
write_csv(x, "compare_filtering.csv")




# ------

x <- cd8_and_tcr_top %>%
    dplyr::left_join(bb_top,
                     by = "tcr_name",
                     suffix = c("_cd8_and_tcrb", "_one_beta")) %>%
    dplyr::relocate(tcr_name,
                    matches("n_clone"),
                    matches("n_cluster"),
                    matches("pct_clone_cluster")) 

write_csv(x, "compare_filtering_strategies.csv")


# Clone tables

x <- cd8_or_tcr %>%
    dplyr::filter(seurat_clusters == 3, ! beta_aa == "NA_NA") %>%
    dplyr::group_by(Sample, beta_aa) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::arrange(desc(n))
write_csv(x, file = "cd8_or_tcr_cluster_3_clones.csv")

y <- one_tcr_beta %>%
    dplyr::filter(seurat_clusters == 1, ! beta_aa == "NA_NA") %>%
    dplyr::group_by(Sample, beta_aa) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::arrange(desc(n))
write_csv(x, file = "one_tcr_beta_cluster_1_clones.csv")

anti_join(x, y, by = c("Sample", "beta_aa")) %>%
    write_csv("in_cd8_or_tcr_not_one_tcr_beta.csv")

anti_join(y, x, by = c("Sample", "beta_aa")) %>%
    write_csv("in_one_tcr_beta_not_cd8_or_tcr.csv")





data_dir <- file.path(workdir, "data/processed")
unfiltered <- file.path(data_dir, "unfiltered/seurat.csv.gz")
md <- read_csv(unfiltered) %>%
    dplyr::mutate(cdr3_aa_beta = CTaa2) %>%
    tidyr::separate(TCR2,
                    into = c("trbv", "trbd", "trbj", "trbc"),
                    sep = "\\.",
                    remove = FALSE) %>%
    dplyr::left_join(selected_clones,
                     by = c("Sample","cdr3_aa_beta", "trbv", "trbj"),
                     relationship = "many-to-one")




# ---------
