# Libraries and resources ----

library("janitor")
library("readxl")
library("tidyverse")

clone_dir <- "data/raw/clones_of_interest/"

neuron_reactive <- file.path(clone_dir, "20250227_List_%s.xlsx")

#--------------------------------------------------------------------------

process_sheet <- function(i, sheet_name, fname){
    new_colnames <- c(patient = "donor",
                      cdr3_aa_beta = "cdr3_aaseq")
    
    readxl::read_excel(fname, sheet = i)  %>%
        janitor::remove_empty(which = c("rows", "cols")) %>%
        janitor::clean_names() %>%
        dplyr::rename(any_of(new_colnames)) %>%
        dplyr::mutate(tcr_category = sheet_name) 
}


process_neuron_rxt <- function(fname_template, markers = c("I","J","K","L")){
    
    all_clones <- lapply(markers, function(nm){
        fname <- sprintf(fname_template, nm)
        sheet_names <- readxl::excel_sheets(fname)
        sheet_names <- gsub(" ", "_", sheet_names)
    
        clones <- lapply(seq_along(sheet_names), function(i){ 
            process_sheet(i, sheet_names[[i]], fname)
        }) %>%
            bind_rows() %>%
            dplyr::mutate(clone_set = nm)
    }) %>% bind_rows()
    
        #dplyr::mutate(Sample = gsub("^HD", "HC", Sample),
       #               tcrb_seq = gsub("(TRB[VJ])0", "\\1", tcrb_seq)) %>%
     
}

#--------------------------------------------------------------------------
md <- read_csv("~/Analyses/Jones_tcell/archive/Feb_2025/data/processed/integrated_seurat.csv.gz")
md <- md %>% dplyr::mutate(patient = gsub("_dis|_5m", "", Sample))
# Add patient to md
