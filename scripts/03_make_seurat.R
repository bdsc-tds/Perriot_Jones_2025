# ----------------------------------------------------------------------------
# Libraries and setup ----

library("hdf5r")
library("rhdf5")
library("scRepertoire")
library("Seurat")
library("SeuratDisk")
library("tidyverse")
library("argparse")

set.seed(7620)

# Command line arguments ----
parser <- ArgumentParser(description = "Create a seurat file")

parser$add_argument('--tcr', '-t', help = 'Input TCR file')
parser$add_argument('--samples', '-s',
                    help = 'Text file with sample names in the first column')
parser$add_argument('--output', '-o', help = 'name for seurat rds output')
parser$add_argument('--cellbender', '-c',
                    help = 'Directory for cellbender input')
parser$add_argument('--min-rna', '-r', dest = "min_rna", 
                    help = 'Minimum total RNA count per cell', default = 200)
parser$add_argument('--clones', dest = "coi", 
                    help = 'Sequence(s) for clones of interest')

clones = "scripts/clones_of_interest.R"
args <- parser$parse_args()

# ----------------------------------------------------------------------------
# Functions ----
# beta_aa ----
beta_aa <- function(md){
    md %>%
        tidyr::separate(CTgene,
                        into = c("TCR1", "TCR2"),
                        sep = "_", remove = FALSE) %>%
        tidyr::separate(CTaa,
                        into = c("CTaa1", "CTaa2"),
                        sep = "_", remove = FALSE) %>%
        dplyr::mutate(across(matches("TCR|CTaa"), ~na_if(.x, "NA")),
                      beta_aa = paste(TCR2, CTaa2, sep = "_")) 
}


# tcr_presence ----
tcr_presence <- function(md){
    md %>%
        dplyr::mutate(
            n_alphas = case_when(is.na(CTaa1) ~ NA,
                                 TRUE ~ stringr::str_count(CTaa1, ";") + 1),
            n_betas = case_when(is.na(CTaa2) ~ NA,
                                TRUE ~ stringr::str_count(CTaa2, ";") + 1),
            tcr_presence = 
                case_when(is.na(CTaa1) & is.na(CTaa2) ~ "no_tcr",
                          n_alphas == 2 & n_betas == 2 ~ "two_chains",
                          n_betas == 2 & is.na(CTaa1) ~ "two_betas_no_alpha",
                          n_betas == 2 & n_alphas == 1 ~ "two_betas_one_alpha",
                          n_alphas == 2 & is.na(CTaa2) ~ "two_alphas_no_beta",
                          n_alphas == 2 & n_betas == 1 ~ "two_alphas_one_beta",
                          n_betas == 1 & is.na(CTaa1) ~ "no_alpha_one_beta",
                          n_alphas == 1 & is.na(CTaa2) ~ "one_alpha_no_beta",
                          n_alphas == 1 & n_betas == 1 ~ "alpha_beta",
                          TRUE ~ "other")) %>%
        return()
}

# CD8 expression ----
cd8_expr <- function(seurat_obj, features = c("CD8A", "CD8B")){
    cd8_expr <- 
        
        seurat_obj[["RNA"]]$counts[features, , drop = FALSE]
    cd8_expressed <- colSums(cd8_expr) > 0
    return(cd8_expressed)
}

# add_tcr_metadata ----
add_tcr_metadata <- function(seurat_obj, combined_tcr){
    # Add TCR information 
    seurat_obj <- combineExpression(combined_tcr, seurat_obj)
    seurat_obj
}

add_beta_aa <- function(seurat_obj){
    # Add beta chain cdr3
    seurat_obj[[]] <- beta_aa(seurat_obj[[]])
    
    # Add vj_aa
    #seurat_obj[[]] <- vj_aa(seurat_obj[[]])
}

# h5_to_seurat ----
h5_to_seurat <- function(nm,
                         exp_dir,
                         min_cells,
                         min_features,
                         keep_rn){
    
    h5_data <- hdf5r::H5File$new(exp_dir, mode = 'r')
    
    decont_m <- Matrix::sparseMatrix(
        i = h5_data[['matrix/indices']][],
        p = h5_data[['matrix/indptr']][],
        x = h5_data[['matrix/data']][],
        dimnames = list(
            h5_data[['matrix/features/name']][],
            h5_data[['matrix/barcodes']][]
        ),
        dims = h5_data[['matrix/shape']][],
        index1 = FALSE
    )
    
    # Add sample to cell names for easier integration of TCR
    colnames(decont_m) <- paste(nm, colnames(decont_m), sep = "_")
    
    samples <- data.frame(Sample = rep(nm, ncol(decont_m)),
                          row.names = colnames(decont_m))
    
    # Create Seurat objects, applying filtering for cells and features
    seurat_obj <- CreateSeuratObject(counts = decont_m[keep_rn, ], 
                                     min.cells = min_cells,
                                     min.features = min_features,
                                     meta.data = samples)
}    


# cellbender_to_seurat ----
# min_cells 10 keep features expressed in at least 10 cells
# min_features 100 require at least 100 genes expressed
# max_percent_mt 10 keep if not more than 10% mitochondrial reads
cellbender_to_seurat <- function(cellbender_dir,
                                 samples,
                                 min_cells = 10,
                                 min_features = 100,
                                 max_percent_mt = 10){
    # Read cellbender filtered data 
    
    exp_dirs <- list.files(cellbender_dir,
                           full.names = TRUE,
                           pattern = paste(samples, collapse = "|"))
    
    names(exp_dirs) <- gsub("_cellbender.*", "", basename(exp_dirs))
    
    # Make sure all samples are found
    stopifnot(all(samples %in% names(exp_dirs)))
    
    exp_dirs <- exp_dirs[samples]
    
    # Get duplicated row names to drop
    # (note that after cellbender, samples do not have the same features)
    rn <- H5File$new(exp_dirs[[1]])[['matrix/features/name']][]
    rhdf5::h5closeAll()
    keep_rn <- ! ( duplicated(rn) | duplicated(rn, fromLast = TRUE) )
    
    # Create Seurat objects from cellbender output
    seurat_objs <- lapply(names(exp_dirs), function(nm){
        print(nm)
        obj <- h5_to_seurat(nm,
                            exp_dirs[[nm]],
                            min_cells,
                            min_features,
                            keep_rn)
        
        obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj,
                                                            pattern = "^MT-")
        
        # Filter high mitochondrial read percentage 
        obj <- subset(obj, subset = percent.mt < max_percent_mt)
        obj
    })
    # Prep for Seurat SketchData -merge and normalize
    seurat_obj <- merge(seurat_objs[[1]], seurat_objs[2:length(seurat_objs)])
    seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
    seurat_obj <- NormalizeData(seurat_obj) # normalise for total count per cell
    seurat_obj
}

#count_tcr_presence ----
count_tcr_presence <- function(md){
    tcrb_exp <- md %>%
        dplyr::group_by(Sample) %>%
        dplyr::mutate(has_tcr = ! tcr_presence == "no_tcr") %>%
        dplyr::summarise(n_sample = n(),
                         cd8b = sum(CD8B),
                         cd8a_or_cd8b = sum(CD8),
                         cd8b_or_tcr = sum(CD8B | has_tcr),
                         cd8b_and_tcr = sum(CD8B & has_tcr),
                         cd8a_or_cd8b_or_tcr = sum(CD8|has_tcr),
                         tcr = sum(has_tcr),
                         tcrb = sum(! is.na(CTaa2)),
                         cd8b_or_tcrb = sum(CD8B | ! is.na(CTaa2)),
                         cd8b_and_tcrb = sum(CD8B & ! is.na(CTaa2))) %>%
        dplyr::mutate(across(.cols = -c("Sample", "n_sample"),
                             .fns = ~ .x/n_sample*100,
                             .names = "{.col}_pct"))
    tcrb_exp
} 


# make_seurat ----
make_seurat <- function(args){
    if (! file.exists(dirname(args$output))){ dir.create(dirname(args$output)) }
    stopifnot(grepl("\\.rds$", args$output))
    
    # Select the samples to use
    samples <- read_delim(args$samples, col_names = FALSE)[[1]]
    
    # Read TCR data
    combined_tcr <- read_rds(args$tcr)
    
    # Read clone of interest
    source(args$coi)
    coi <- get_coi()
    
    # Check that sample names are the same as in the config file
    stopifnot(length(intersect(names(combined_tcr), samples)) ==
                 length(samples))
    
    # Create Seurat object with minimal filtering
    seurat_obj <- cellbender_to_seurat(args$cellbender, samples)
    
    # Add disease to metadata 
    seurat_obj$condition <- gsub("[0-9]+$", "", seurat_obj$Sample)

    # Add TCR information
    seurat_obj <- add_tcr_metadata(seurat_obj, combined_tcr)
    seurat_obj[[]] <- add_beta_aa(seurat_obj)
    
    # Add CD8 expression status to metadata
    seurat_obj$CD8 <- cd8_expr(seurat_obj)
    seurat_obj$CD8B <- cd8_expr(seurat_obj, features = "CD8B")
    
    # Classify cells by tcr presence
    seurat_obj[[]] <- tcr_presence(seurat_obj[[]])
    
    # Filter for minimum RNA expressed
    #seurat_obj <- subset(seurat_obj, nCount_RNA >= args$min_rna)
    
     # Add clone of interest to metadata (must be after tcr metadata added)
    seurat_obj[[]] <- seurat_obj[[]] %>%
        dplyr::mutate(coi = case_when(beta_aa == coi$ri01_coi ~ "Ri01", 
                      beta_aa == coi$ri02_coi ~ "Ri02",
                      TRUE ~ "no"))
    
    seurat_obj$is_coi = ifelse(seurat_obj$coi == "no", "non-coi", "coi")
    
    # Save seurat object and metadata
    write_rds(seurat_obj, args$output)
    write_csv(seurat_obj[[]], gsub("rds$", "csv.gz", args$output))
    
    tcrb_exp <- count_tcr_presence(seurat_obj[[]])
    write_csv(tcrb_exp, gsub(".rds$", "_tcr_presence.csv", args$output))
    
    return()
}

# ----------------------------------------------------------------------------   
# Run ----

make_seurat(args)
