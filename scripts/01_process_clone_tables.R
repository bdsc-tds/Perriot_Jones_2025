# Libraries and resources ----

library("janitor")
library("readxl")
library("tidyverse")

clone_dir <- "data/raw/clones_of_interest/"

neuron_reactive <- file.path(clone_dir, "20250227_List_K.xlsx")
selected_clones <- file.path(clone_dir, "20250228_List_%s.xlsx")
selected_clones <- sprintf(selected_clones, "N")

clones_w_alpha <- "data/raw/20250114_TCR_list_Helen.xlsx"
processed_clones_out <- "data/processed/merged_clone_sets.csv" 

#--------------------------------------------------------------------------

# Patch alpha chain ----
alpha_patch <- function(){
    tibble::tribble(
        ~patient, ~cdr3_aa_alpha, ~cdr3_aa_beta, ~trav, ~traj, ~trac,
        "Ri01", "CAASSTGGFKTIF", "CASSPGLGNYEQYF", "TRAV13-1", "TRAJ9", "TRAC",
        "Ri02", "CAFMNPRRDDKIIF", "CASSLERGLGGYTF", "TRAV38-1", "TRAJ30", "TRAC")
}
    
patch_alpha <- function(clones, patch_cols = NULL){
    
    patch <- alpha_patch()
    if (is.null(patch_cols)){
        patch_cols <- intersect(colnames(patch), colnames(clones))
    }
    # Check that patch_cols columns are sufficient to uniquely identify clones
    check_id_cols <- clones %>%
                     dplyr::group_by(across(all_of(patch_cols))) %>%
                     dplyr::summarise(n = n()) %>%
                     dplyr::pull(n) %>%
                     max()
    if (check_id_cols > 1) {
        stop("id columns don't uniquely identify clones")
    } 
    
    clones <- clones %>%
        dplyr::rows_patch(patch %>% dplyr::select(all_of(patch_cols))) 
    
    return(clones)
}

# Process single sheet using read_func function that only accepts fname ----
process_sheet <- function(read_func, sheet_name, fname){
    new_colnames <- c(patient = "donor",
                      cdr3_aa_beta = "cdr3_aaseq")
    
    read_func(fname) %>%    
        janitor::remove_empty(which = c("rows", "cols")) %>%
        janitor::clean_names() %>%
        dplyr::rename(any_of(new_colnames)) %>%
        dplyr::mutate(tcr_category = sheet_name) 
}

# Process all sheets / csv files in a directory ----
process_all_sheets <- function(fname){
    # If it's directory of csvs:
    if (dir.exists(fname)){
        fnames <- list.files(fname, pattern = ".csv", full.names = TRUE)
        sheet_names <- gsub("-Table.*", "", basename(fnames))
        clones <- lapply(seq_along(sheet_names), function(i){
            process_sheet(partial(read_csv,
                                  locale = readr::locale(encoding = "UTF-8")),
                          sheet_names[i], fnames[i]) 
        }) 
    } else {
        
        # If it's an Excel sheet:
        sheet_names <- readxl::excel_sheets(fname)
        sheet_names <- gsub(" ", "_", sheet_names)
            
        clones <- lapply(seq_along(sheet_names), function(i){ 
            process_sheet(partial(read_excel, sheet = i),
                          sheet_names[[i]], fname)
        }) 
    }
    return(clones %>% dplyr::bind_rows())
}

# Reformat tcr name ----
reformat_tcr_name <- function(df){
    df %>% 
        dplyr::mutate(tcr_name = gsub('[\"\u00A0]', "", tcr_name),
                      tcr_name = gsub('\\s+', " ", tcr_name),
                      tcr_name = gsub("-TCR", "- TCR", tcr_name),
                      tcr_name = gsub("((HD|Ri)[0-9]+)- TCR", "\\1 - TCR", tcr_name),
                      
                      tcr_name_alpha = gsub("\\s", "", tcr_name),
                      # Consolidate (1) and (2) into single entries
                      tcr_name_alpha = gsub("\\(.*", "", tcr_name_alpha),
                      tcr_name_alpha = gsub("-", "_", tcr_name_alpha),
                      tcr_name_alpha = gsub("TCR", "TCR_", tcr_name_alpha))
}

# Process clone table containing alpha chains ---- 
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

# Process neuron reactive clones ----
process_neuron_rxt <- function(reactive_clones,
                               selected_clones,
                               #alpha_clones,
                               outfname){
    # Only list K is needed:
    # List I = all neuron reactive
    # List J = all Ri
    # List K = all clones
    # List L = all HD
    
    # List M = all clones in cluster 2 with at least 5 reads
    # List N = selected clones 
    
    # Selected clones for differential expression analyses ----
    selected_clones <- process_all_sheets(selected_clones) %>%
        dplyr::mutate(selected_clone = TRUE,
                      Sample = case_when(patient == "Ri01" ~ "Ri01_dis",
                                      TRUE ~ patient)) %>%
        dplyr::select(-tcr_category) %>%
        reformat_tcr_name() 

    clones <- process_all_sheets(reactive_clones) %>%
        dplyr::mutate(tcr_category = gsub("_.*", "", tcr_category),
                      reactive = 
                          case_when(tcr_category == "Neuron-reactive" ~ TRUE,
                                    tcr_category == "Non-reactive" ~ FALSE,
                                    TRUE ~ NA)) %>%
        dplyr::select(-tcr_category) %>%
        reformat_tcr_name() %>%
        
        # Manually checked that this doesn't result in extra rows that
        # should be merged
        full_join(selected_clones) 
    
   # Merge clones sharing beta chain
    clones <-
        clones %>%
        dplyr::group_by(patient, tcr_name_alpha) %>%
        # Classify as reactive if any clone is reactive
        dplyr::mutate(reactive = any(reactive == TRUE)) %>%
        dplyr::select(-tcr_name) %>%
        dplyr::rename(tcr_name = tcr_name_alpha) %>%
        dplyr::ungroup() %>%
        # Fix extra 0 in trbj
        dplyr::mutate(trbj = gsub("TRBJ0", "TRBJ", trbj)) %>%
        unique()
    
    # Set N contains clones from a variety of clusters, mostly cluster 2
    write_csv(clones, file = outfname)
    invisible(clones)
}

#--------------------------------------------------------------------------

process_neuron_rxt(reactive_clones = neuron_reactive,
                   selected_clones = selected_clones,
                   outfname = processed_clones_out)

