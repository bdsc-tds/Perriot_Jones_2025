# Libraries and resources ----

library("janitor")
library("readxl")
library("tidyverse")

clone_dir <- "data/raw/clones_of_interest/"

neuron_reactive <- file.path(clone_dir, "20250227_List_K.xlsx")
selected_clones <- file.path(clone_dir, "20250228_List_%s.xlsx")

clones_w_alpha <- "data/raw/20250114_TCR_list_Helen.xlsx"



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

process_all_sheets <- function(nm, fname){
    sheet_names <- readxl::excel_sheets(fname)
    sheet_names <- gsub(" ", "_", sheet_names)
        
    clones <- lapply(seq_along(sheet_names), function(i){ 
        process_sheet(i, sheet_names[[i]], fname)
    }) %>%
        bind_rows() %>%
        dplyr::mutate(clone_set = nm)
    
    return(clones)
}


process_alpha_clones <- function(fname){
    readxl::read_excel(fname, sheet = 1)  %>%
        janitor::remove_empty(which = c("rows", "cols")) %>%
        janitor::clean_names() %>%
        dplyr::mutate(sample = gsub("_.*", "", name)) %>%
        dplyr::relocate(sample) %>%
        
        # Rename the CDR3 columns 
        dplyr::rename("cdr3_aa_alpha" = cdr3_aaseq_3,
                      "cdr3_aa_beta" = cdr3_aaseq_6,
                      "patient" = sample,
                      "tcr_name_alpha" = name) %>%
        
        dplyr::mutate(alpha_vj_aa = paste(trav, traj, cdr3_aa_alpha, sep = "."),
                      beta_vj_aa = paste(trbv, trbj, cdr3_aa_beta, sep = "."),
                      vj_aa = paste(alpha_vj_aa, beta_vj_aa, sep = "_")) 
}


process_neuron_rxt <- function(reactive_template,
                               selected_clones,
                               markers = c("K","M","N")){
    # List I = all neuron reactive
    # List J = all Ri
    # List K = all clones
    # List L = all HD
    
    # List with alpha chains
    alpha_clones <- process_alpha_clones(clones_w_alpha)
    
    
    all_clones <- lapply(markers, function(nm){
        fname <- sprintf(fname_template, nm)
        sheet_names <- readxl::excel_sheets(fname)
        sheet_names <- gsub(" ", "_", sheet_names)
    
        clones <- lapply(seq_along(sheet_names), function(i){ 
            process_sheet(i, sheet_names[[i]], fname)
        }) %>%
            bind_rows() %>%
            dplyr::mutate(clone_set = nm)
    }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(tcr_category = gsub("_.*", "", tcr_category))
       
    clones_wide <- all_clones %>%
        tidyr::pivot_wider(id_cols = c("patient",
                                       "tcr_name",
                                       "trbv",
                                       "trbj",
                                       "cdr3_aa_beta"),
                           names_from = clone_set,
                           values_from = tcr_category) 
        
    
    
    
     
    x <- all_clones %>%
        left_join(alpha_clones, relationship = "many-to-many")
    
    
    
    all_clones %>%
        dplyr::group_by(across(all_of(c("patient", "tcr_name","trbv", "trbj")))) %>%
        dplyr::filter(n_distinct(tcr_category) > 1) 
    
     
        #dplyr::mutate(Sample = gsub("^HD", "HC", Sample),
       #               tcrb_seq = gsub("(TRB[VJ])0", "\\1", tcrb_seq)) %>%
     
}

#--------------------------------------------------------------------------
md <- read_csv("~/Analyses/Jones_tcell/archive/Feb_2025/data/processed/integrated_seurat.csv.gz")
md <- md %>%
    dplyr::mutate(patient = gsub("_dis|_5m", "", Sample)) %>%
    dplyr::select(patient, matches("CT|TCR"), beta_aa, vj_aa) %>%
    unique()



# The original clones with alpha that aren't found
md %>% filter(CTaa2 == "CASSPGLGNYEQYF", patient == "Ri01")
# TRAV13-1.TRAJ9.CAASSTGGFKTIF_TRBV28.TRBJ2-7.CASSPGLGNYEQYF
# NA.NA.NA_TRBV28.TRBJ2-7.CASSPGLGNYEQYF

# one clone missing trbc
#x <- md %>%
#    dplyr::filter(CTaa2 == "CASSLERGLGGYTF", patient == "Ri02") %>%
#    data.frame() 

#TRAV38-1.TRAJ30.CAFMNPRRDDKIIF_TRBV13.TRBJ1-2.CASSLERGLGGYTF
#NA.NA.NA_TRBV13.TRBJ1-2.CASSLERGLGGYTF



