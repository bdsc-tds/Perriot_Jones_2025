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

# Make shortened ID column vj_aa ----
vj_aa <- function(md){
    temp <- md %>%
        dplyr::select(TCR1, TCR2, CTaa1, CTaa2) %>%
        dplyr::mutate(across(c(TCR1, TCR2, CTaa1, CTaa2), ~strsplit(.x, ";")),
                      id = 1:n()) %>%
        tidyr::unnest(c(TCR1, TCR2, , CTaa1, CTaa2)) %>%
        tidyr::separate_wider_delim(cols = c(TCR1, TCR2),
                                    delim = '.', names_sep = "_") %>%
        dplyr::mutate(beta_vj_aa = paste(TCR2_1, TCR2_3, CTaa2, sep = "."),
                      alpha_vj_aa = paste(TCR1_1, TCR1_2, CTaa1, sep = "."),
                      vj_aa = paste(alpha_vj_aa, beta_vj_aa, sep = "_")) %>%
        dplyr::select(id, vj_aa) %>% 
        dplyr::group_by(id) %>%
        dplyr::mutate(row_n = 1:n()) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = "row_n",
                           names_prefix = "tcr_",
                           values_from = "vj_aa") %>%
        # If there are multiple TCRs
        tidyr::unite(col = vj_aa, tcr_1, tcr_2, sep = ";") %>%
        dplyr::mutate(vj_aa = gsub(";NA$", "", vj_aa),
                      vj_aa = na_if(vj_aa, "NA.NA.NA_NA.NA.NA")) %>%
        dplyr::select(-id)
    md <- md  %>% dplyr::bind_cols(temp)
    return(md)
}

# tcr_presence ----
tcr_presence <- function(md){
    md %>%
        dplyr::mutate(
            two_alphas = grepl(";", CTaa1),
            two_betas = grepl(";", CTaa2),
            two_chains = two_alphas & two_betas,
            tcr_presence = 
                case_when(is.na(CTaa1) & is.na(CTaa2) ~ "no_alpha_no_beta",
                          two_chains ~ "two_chains",
                          two_betas & is.na(CTaa1) ~ "two_betas_no_alpha",
                          two_betas & !is.na(CTaa1) ~ "two_betas_one_alpha",
                          two_alphas & is.na(CTaa2) ~ "two_alphas_no_beta",
                          two_alphas & !is.na(CTaa2) ~ "two_alphas_one_beta",
                          is.na(CTaa1) ~ "no_alpha",
                          is.na(CTaa2) ~ "no_beta",
                          ! is.na(CTaa1) & ! is.na(CTaa2) ~ "alpha_beta",
                          TRUE ~ "unknown")) %>%
        dplyr::select(-two_alphas, -two_betas) %>%
        return()
}

# CD8 expression ----
cd8_expr <- function(seurat_obj){
    cd8_expr <- subset(seurat_obj, features = c("CD8A", "CD8B"))
    cd8_expressed <- colSums(cd8_expr@assays$RNA$counts) > 0
    return(cd8_expressed)
}

# add_tcr_metadata ----
add_tcr_metadata <- function(seurat_obj, combined_tcr){
    # Add TCR information 
    seurat_obj <- combineExpression(combined_tcr, seurat_obj)
    
    # Add beta chain cdr3
    seurat_obj[[]] <- beta_aa(seurat_obj[[]])
    
    # Add vj_aa
    seurat_obj[[]] <- vj_aa(seurat_obj[[]])
    
    # Add disease to metadata 
    seurat_obj$condition <- gsub("[0-9]+$", "", seurat_obj$Sample)
    seurat_obj
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
#
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

# make_seurat ----
make_seurat <- function(args){
    if (! file.exists(dirname(args$output))){ dir.create(dirname(args$output)) }
    stopifnot(grepl("\\.rds$", args$output))
    
    # Select the samples to use
    samples <- read_delim(args$samples, col_names = FALSE)[[1]]
    
    # Read TCR data
    combined_tcr <- read_rds(args$tcr)
    
    # Check that sample names are the same as in the config file
    stopifnot( length(intersect(names(combined_tcr), samples)) ==
                   length(samples))
    
    # Create Seurat object with minimal filtering
    seurat_obj <- cellbender_to_seurat(args$cellbender, samples)
    
    # Filter for minimum RNA expressed
    seurat_obj <- subset(seurat_obj,
                         cells = Cells(seurat_obj)[seurat_obj$nCount_RNA >=
                                                       args$min_rna])
    
    # Add TCR information
    seurat_obj <- add_tcr_metadata(seurat_obj, combined_tcr)
    
    # Add CD8 expression status to metadata
    seurat_obj$CD8 <- cd8_expr(seurat_obj)
    
    # Classify cells by tcr presence
    seurat_obj[[]] <- tcr_presence(seurat_obj[[]])
    
    # Save seurat object and metadata
    write_rds(seurat_obj, args$output)
    write_csv(seurat_obj[[]], gsub("rds$", "csv.gz", args$output))
    
    return()
}

# ----------------------------------------------------------------------------   
# Run ----

make_seurat(args)
