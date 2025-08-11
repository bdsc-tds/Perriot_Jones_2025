#!/bin/bash

module load singularityce

# ENVIRONMENT VARIABLES
WORKDIR=${WORKDIR}
CONTAINER_DIR=${CONTAINER_DIR}
SCRATCH_DIR=${SCRATCH_DIR}

export SINGULARITY_BINDPATH="/users,/scratch,/work,/data"

singularity shell ${CONTAINER_DIR}tcr_seq.sif

# Output directories
FILTER_NM="one_tcr_beta"

DATA="${WORKDIR}data/processed/"
FILTER_DATA="${WORKDIR}data/processed/${FILTER_NM}/"

RESULTS="${WORKDIR}results/${FILTER_NM}"
SCRIPTS="${WORKDIR}scripts/"

# Other inputs
TCR="${WORKDIR}data/processed/combined_tcr.rds"
TCR_CATS="${DATA}tcr_sequences.csv"
CB_INPUT="${WORKDIR}data/cellbender_input.txt"
CB_DIR="${WORKDIR}data/processed/cellbender/filtered"
INTEGRATION_METHOD="rpca"

#--------------------------------------------------------------------------
# Combine cellbender samples into Seurat object
Rscript ${SCRIPTS}02_make_seurat.R \
--tcr ${TCR} --samples ${CB_INPUT} \
--output ${DATA}unfiltered/seurat.rds \
--cellbender ${CB_DIR} \
--clones ${SCRIPTS}clones_of_interest.R

# Apply filters
Rscript ${SCRIPTS}03_filter_seurat.R \
--filter ${FILTER_NM} \
--seurat ${DATA_OUT}"unfiltered/seurat.rds"\
--output ${FILTER_DATA} 

# Run integration
Rscript ${SCRIPTS}04_run_integration.R \
--input ${FILTER_DATA_OUT}filtered_seurat.rds \
--output ${FILTER_DATA_OUT} \
--figures ${RESULTS} \
--integration ${INTEGRATION_METHOD}

# fig_4_bc_supp_fig_1a_2a (originally global_analyses_umaps)
Rscript ${SCRIPTS}fig_4_bc_supp_fig_1a_2a.R \
--seurat ${FILTER_DATA}integrated_seurat.rds \
--results ${RESULTS}

# supp_fig_1b 
Rscript ${SCRIPTS}supp_fig_1b.R \
--seurat ${FILTER_DATA}integrated_seurat.rds \
--results ${RESULTS}


# supp_fig_2b_2c.R 
Rscript ${SCRIPTS}supp_fig_2b_2c.R \
--metadata ${FILTER_DATA}integrated_seurat.csv.gz \
--results ${RESULTS}

# figure_4f.R 
Rscript ${SCRIPTS}figure_4f.R \
--seurat ${FILTER_DATA}integrated_seurat.rds \
--results ${RESULTS} \
--workdir ${WORKDIR}

# supp_table_6_diff_expr_clusters.R 
Rscript ${SCRIPTS}supp_table_6_diff_expr_clusters.R \
--seurat ${FILTER_DATA}integrated_seurat.rds \
--results ${RESULTS} 

