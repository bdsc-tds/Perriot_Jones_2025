# Libraries and resources ----

library("janitor")
library("readxl")
library("tidyverse")

clone_dir <- "data/raw/clones_of_interest/"

neuron_reactive <- file.path(clone_dir, "20250227_List_K.xlsx")
selected_clones <- file.path(clone_dir, "20250228_List_%s.xlsx")
selected_clones <- sprintf(selected_clones, "N")

clones_w_alpha <- "data/raw/20250114_TCR_list_Helen.xlsx"
processed_clones_out <- "data/processed/clone_sets.csv" 

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
    
    # The original clones with alpha that aren't found
    #md %>%
    #    dplyr::filter(CTaa2 == "CASSLERGLGGYTF", patient == "Ri02") %>%
    #    tidyr::separate(TCR1, sep = "\\.", into = c("trav", "traj", "trac")) %>%
    #    dplyr::rename("cdr3_aa_alpha" = CTaa1,
    #                  "cdr3_aa_beta" = CTaa2) %>%
    #    dplyr::filter(! is.na(trav)) %>%
    #    dplyr::select(patient, matches("^cdr3_aa"), matches("^tra")) %>%
    #    unique()
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

reformat_tcr_name <- function(df){
    df %>% 
        dplyr::mutate(tcr_name = gsub('[\"\u00A0]', "", tcr_name),
                      tcr_name = gsub('\\s+', " ", tcr_name),
                      tcr_name = gsub("-TCR", "- TCR", tcr_name),
                      tcr_name = gsub("((HD|Ri)[0-9]+)- TCR", "\\1 - TCR", tcr_name),
                      
                      tcr_name_alpha = gsub("\\s", "", tcr_name),
                      tcr_name_alpha = gsub("\\(", " \\(", tcr_name_alpha),
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
                               alpha_clones,
                               outfname){
    # Only list K is needed:
    # List I = all neuron reactive
    # List J = all Ri
    # List K = all clones
    # List L = all HD
    
    # List M = all clones in cluster 2 with at least 5 reads
    # List N = selected clones 
    
    # List with alpha chains
    alpha_clones <- process_alpha_clones(alpha_clones)
    
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
        
        
    # Name is needed to add the alpha 
    # Patch name to match name used reactivity 
    no_tcr_name <- clones %>%
        dplyr::filter(is.na(tcr_name))
        
    clones <- clones %>%    
        dplyr::filter(! is.na(tcr_name)) %>%

        # Not all alpha clones are in the reactive clones
        # Two reactive clones are missing from the alpha clones, these will be
        # patched
        dplyr::left_join(alpha_clones, relationship = "one-to-one") %>%
        
        # Only two reactive clones are not in the alpha clones, 
        # these will be patched 
        patch_alpha() %>%
        
        # Add in the unmatched tcrs again
        dplyr::bind_rows(no_tcr_name)
  
    # Set N contains clones from a variety of clusters, mostly cluster 2
    write_csv(clones, file = outfname)
    invisible(clones)
}

# Scratch ----
scratch(){
    # Code used for checking, not run in pipeline
    
    # Check that all of the selected clones appear in the meta data ----
    found_by_sample <- dplyr::semi_join(selected_clones, md,
                                        by = c("Sample", "cdr3_aa_beta"))
    
    found_by_patient <- dplyr::semi_join(selected_clones, md,
                                         by = c("patient", "cdr3_aa_beta"))
    
    missing <- selected_clones %>%
        dplyr::anti_join(found_by_sample) %>%
        dplyr::anti_join(found_by_patient)
    
    # Check that selected cluster two is in cluster 2 ----
    selected_cl2 <- selected_clones %>%
        dplyr::filter(tcr_category =="Cluster_2_Ri") 
    
    selected_cl2 %>%
        dplyr::anti_join(md %>% filter(seurat_clusters == "2"),
                         by = c("Sample", "cdr3_aa_beta"))
    
    # One clone isn't in Cl 2 in this metadata 
 #   TRBV4-2 CASAPQGAPGNTIYF TRBJ1-3 Cluster_2_Ri  Ri05 

    # Check if selected cluster two is everything from cluster 2 ----
    md_c2_ri <- md %>%
        dplyr::filter(seurat_clusters == "2", grepl("Ri", Sample),
                      ! grepl("Ri01_5m", Sample))
              
    missing_from_cl2 <- md_c2_ri %>% dplyr::anti_join(selected_cl2, 
                                  by = c("Sample", "cdr3_aa_beta")) %>%
        group_by(Sample, cdr3_aa_beta) %>%
        summarise(n = n()) %>% arrange(desc(n))
    
    # All missing from cl2 have at most 4 reads, were they filtered?
    
    # Set N - everything in set N is in the metadata
    set_n %>% anti_join(md, by = c("patient", "cdr3_aa_beta"))
    
    # Set N also needs the alpha chains for distinguishing
    set_n %>% left_join(md, by = c("patient", "cdr3_aa_beta"))
    
    # Check if set_N is in the reactive clones set - two are missing ----
    #Ri02   TRBV13       CASSLAGYNSPLHF TRBJ1-6 Ri02 - TCR 1  
    #HD7 TRBV20-1 CSARAPGTGGWLRGYQPQHF TRBJ1-5  HD7 -  TCR 2  
 
    # Two are missing from the alpha clones (same two as reactive)
    n_missing <- set_n %>% anti_join(alpha_clones, 
                                     by = c("patient", "cdr3_aa_beta"))
    # Ri01 TRBV28 CASSPGLGNYEQYF TRBJ02-7 "old TCR 3"  
    # Ri02 TRBV13 CASSLERGLGGYTF TRBJ01-2 "TCR C"
    
    # Check if they have a single alpha chain in the MD
    md %>%
        dplyr::semi_join(n_missing, by = c("patient", "cdr3_aa_beta")) %>%
        dplyr::select(-seurat_clusters, -CTstrict, -CTnt, -Sample) %>%
        unique() %>%
        data.frame() 

    # Check that the selected clones are in the metadata
    
    found_by_sample <- dplyr::semi_join(selected_clones, md,
                                        by = c("Sample", "cdr3_aa_beta"))
    
    found_by_patient <- dplyr::semi_join(selected_clones, md,
                                         by = c("patient", "cdr3_aa_beta"))
    
    missing <- selected_clones %>%
        dplyr::anti_join(found_by_sample) %>%
        dplyr::anti_join(found_by_patient)
    
    
    dplyr::anti_join(selected_clones, reactive_clones,
                     by = c("patient", "cdr3_aa_beta"))
    
    # Two of the selected clones with tcr names are not in the reactivity table
    #Ri02   TRBV13       CASSLAGYNSPLHF TRBJ1-6 Selected_clones_Ri        
    #HD7 TRBV20-1 CSARAPGTGGWLRGYQPQHF TRBJ1-5 Selected_clones_HD 
    
    # Safe to assume Ri01 refers to dis? 
    
    # Checking whether selected clones should be patched by reactive_clones
    # Unambiguous clones from reactive clones
    reactive_patch <- reactive_clones %>%
        dplyr::group_by(patient, cdr3_aa_beta) %>%
        dplyr::filter(n() == 1)

    selected_clones %>% 
        dplyr::filter(! is.na(tcr_name)) %>%
        semi_join(reactive_patch, by = c("patient","cdr3_aa_beta")) %>%
        anti_join(reactive_patch, by = c("patient","cdr3_aa_beta",
                                         "tcr_name_alpha"))    
    
    # Check if alpha chains are unambiguous for missing TCR names
    x <- selected_clones %>%
        dplyr::filter(is.na(tcr_name)) 
    
    alpha_clones %>%
        semi_join(x, by = c("patient", "cdr3_aa_beta")) %>%
        group_by(patient, cdr3_aa_beta) %>%
        filter(n() > 1)
    # Need the names to match the alpha chains.
    
    # (cut from processing function, for processing set M and N together)
    # Selected clones for differential expression analyses ----
    selected_clones <- lapply(sets, function(nm){
        fname <- sprintf(selected_clones, nm)
        process_all_sheets(fname) %>%
            dplyr::mutate(clone_set = nm)
    }) %>%
        dplyr::bind_rows() %>%
        # Set M contains Sample, not patient 
        dplyr::mutate(Sample = case_when(clone_set == "M"~patient,
                                         TRUE ~ NA_character_),
                      patient = case_when(clone_set == "N" ~patient,
                                          clone_set == "M" ~
                                              gsub("_.*", "", Sample))) %>%
        # Fill in patient (assuming there are no clones from Ri01_5m) 
        dplyr::group_by(patient) %>%
        tidyr::fill(Sample, .direction = "updown") %>%
        dplyr::ungroup() %>%
        # Set name not needed anymore as it is coded in tcr_category
        # Cluster 2 not ne
        dplyr::select(-clone_set) %>%
        
        dplyr::mutate(selected_clone =
                          case_when(grepl("Selected_clones", tcr_category ~ TRUE,
                                          TRUE ~ FALSE))) %>%
        
        # Checks on joined table
        clones %>% group_by(tcr_name) %>% filter(n_distinct(cdr3_aa_beta) > 1)
    
       # All reactive clones are in selected clones?
}


#--------------------------------------------------------------------------

process_neuron_rxt(reactive_clones = neuron_reactive,
                   selected_clones = selected_clones,
                   alpha_clones = clones_w_alpha,
                   outfname = processed_clones_out)

#--------------------------------------------------------------------------

reactive_clones <- neuron_reactive


# metadata for checking clones
#md <- read_csv("~/Analyses/Jones_tcell/archive/Feb_2025/data/processed/integrated_seurat.csv.gz")
#md <- md %>%
#    dplyr::mutate(patient = gsub("_dis|_5m", "", Sample)) %>%
#    dplyr::select(Sample, patient, matches("CT|TCR"),
#                   beta_aa, vj_aa, seurat_clusters) %>%
#    unique() %>%
#    dplyr::rename(cdr3_aa_beta = CTaa2)

#

