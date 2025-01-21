# ----------------------------------------------------------------------------
# Libraries and setup ----

library("Seurat")
library("SeuratDisk")
library("tidyverse")
library("argparse")

parser <- ArgumentParser(description = "Create a seurat file")

parser$add_argument('--seruat', '-s', help = 'Input Seurat, as rds')
parser$add_argument('--annotation', '-a', help = 'Table containing cells to annotate')
args <- parser$parse_args()


# ----------------------------------------------------------------------------
# update_seurat ----
update_seurat <- function(seurat_obj, new_md, md_anno, filename){
    new_md <- new_md %>% dplyr::select(-Cell) %>% as.data.frame()
    seurat_obj[[]] <- new_md 
    write_rds(seurat_obj, file = filename)
    write_csv(seurat_obj[[]], file = gsub(".rds$", ".csv.gz", filename))
    write_csv(md_anno, file = gsub(".rds$", "_md_wide.csv.gz", filename))
}

# join_md_tcr ----
join_md_tcr <- function(md, tcr){
    md %>%
        # Only keep rows that match tcr_sequences
        dplyr::inner_join(tcr) %>%
        dplyr::group_by(Cell) %>%
        dplyr::rename(old_name = tcr_name) %>%
        dplyr::mutate(all_names = paste(sort(unique(old_name)), collapse = ", "),
                      
                      # Recode the name
                      tcr_name = case_when(tcr_category == old_cat ~ old_name,
                                           TRUE ~ NA_character_)) %>%
        dplyr::ungroup() 
}

# format_tcr_seqs ----
format_tcr_seqs <- function(tcr_seqs){
    # vj_aa is removed from tcr_seqs so it doesn't contribute to matching
    tcr_seqs %>%
        dplyr::select(sample, 
                      tcr_name, 
                      matches("^cdr3|^tr"))
}

# make_long_md ----
make_long_md <- function(md){
    md %>%
        dplyr::select(matches("[Ss]ample"),
                      matches("Cell|TCR|CTaa[12]|vj_aa")) %>%
        # Both TCRs must be present
        dplyr::filter(! is.na(TCR1), ! is.na(TCR2)) %>%
        dplyr::mutate(across(c(TCR1, TCR2, CTaa1, CTaa2),
                             ~strsplit(.x, ";"))) %>%
        tidyr::unnest(c(TCR1, TCR2, CTaa1, CTaa2)) %>%
        tidyr::separate_wider_delim(cols = c(TCR1, TCR2),
                                    delim = '.', names_sep = "_") %>%
        dplyr::rename(trav = TCR1_1,
                      traj = TCR1_2,
                      trbv = TCR2_1, 
                      trbj = TCR2_3,
                      cdr3_aa_alpha = CTaa1,
                      cdr3_aa_beta = CTaa2) 
}


# classify_tcrs ----
classify_tcrs_md <- function(tcr, md){
    
    tcr_seqs <- read_csv(tcr)
    
    md <- as_tibble(md, rownames = "Cell")
    
    # Subset to just the variables required for joining
    tcr_sub <- format_tcr_seqs(tcr_seqs)
    
    # If md is already annotated, remove the annotations
    md <- md %>%
        select(! any_of(c("tcr_category",
                          "tcr_category_label",
                          "tcr_name",
                          "all_cats",
                          "all_names")))
    
    # Make long format metadata table and join with tcr sequences
    md_long <- make_long_md(md)
    md_anno <- join_md_tcr(md_long, tcr_sub, tcr_name_to_label)
    new_md <- make_wide_md(md, md_anno)
    
    return(list(md_long = md_long,
                md_anno = md_anno,
                new_md = new_md))
}


annotate_tcrs <- function(tcr, seurat_fname){
    seurat_obj <- read_rds(seurat_fname)
    
    res <- classify_tcrs_md(tcr, seurat_obj)
    
    # Add the new metadata table and overwrite the seurat object
    update_seurat(seurat_obj, res$new_md, res$md_anno, seurat_fname)
}

# ----------------------------------------------------------------------------

filter_seurat(args)
    