
workdir='/data/PRTNR/CHUV/DIR/rgottar1/encephalitis_tcr/Perriot_Jones/'
if (! file.exists(workdir)) { workdir <- "." }

tcr <- file.path(workdir, "data/processed/combined_tcr.rds")
cb_input <- file.path(workdir, "data/cellbender_input.txt")
seurat_out <- file.path(workdir, "data/processed/unfiltered/seurat.rds")
cb_dir <- file.path(workdir, "data/processed/cellbender/filtered")
cd8_and_tcr_data <- file.path(workdir, "data/processed/cd8_and_tcrb")
cd8_and_tcr_res <- file.path(workdir, "results/cd8_and_tcrb")
cd8_and_tcr_seurat <- file.path(cd8_and_tcr_data, "filtered_seurat.rds")
cd8_and_tcr_integrated <- file.path(cd8_and_tcr_data, "integrated_seurat.rds")

# 02_make_seurat.R
args <- list(tcr=tcr,
             samples=cb_input,
             output=seurat_out,
             cellbender=cb_dir,
             min_rna=200)

# 03_filter_seurat.R
args <- list(filter="cd8_and_tcr",
             seurat=seurat_out,
             output=cd8_and_tcr_data)

# 04_run_integration.R
args <- list(input=cd8_and_tcr_seurat,
             output=cd8_and_tcr_data,
             figures=cd8_and_tcr_res, # CHECK HERE 
             integration="rpca",
             n_sketch=6000,
             cluster_res=0.5)

#05_barplots.R
args <- list(metadata = file.path(cd8_and_tcr_data, "integrated_seurat.csv.gz"),
             figures = file.path(cd8_and_tcr_res, "cluster_barplots"))

#05_tables.R
args <- list(metadata = file.path(cd8_and_tcr_data, "integrated_seurat.csv.gz"),
             output = file.path(cd8_and_tcr_res, "tables"))

#06_seurat_plots.R
args <- list(metadata = file.path(cd8_and_tcr_data, "integrated_seurat.rds"),
             figures = cd8_and_tcr_res) # FIGURES HERE IS BASEDIR




##08_custom_umaps.R
#args <- list(metadata = file.path(cd8_and_tcr_data, "integrated_seurat.csv.gz"),
#             figures = file.path(cd8_and_tcr_res, "umap"))

## 09_differential_expression.R
#args <- list(input=cd8_and_tcr_integrated,
#             output=file.path(cd8_and_tcr_res, "differential_expression"))
