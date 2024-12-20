#!/bin/bash

WORKDIR='/data/PRTNR/CHUV/DIR/rgottar1/encephalitis_tcr/Perriot_Jones/'

# Output directories
FILTER_NM='cd8_and_tcrb'
DATA_OUT='${WORKDIR}data/processsed/'
RESULTS='${WORKDIR}results/${FILTER_NM}'
FILTER_DATA_OUT='${WORKDIR}data/processsed/${FILTER_NM}/'
SCRIPTS='${WORKDIR}scripts/'
# Other inputs
TCR='${WORKDIR}data/processed/combined_tcr.rds'
CB_INPUT='${WORKDIR}data/cellbender_input.txt'
CB_DIR='${WORKDIR}data/processed/cellbender'
INTEGRATION_METHOD='rpca'

#--------------------------------------------------------------------------
# Combine cellbender samples into Seurat object
Rscript ${SCRIPTS}02_make_seurat.R \
--tcr ${TCR} --samples ${CB_INPUT} \
--output ${DATA_OUT} \
--cellbender ${CB_DIR} 

# Apply filters
Rscript ${SCRIPTS}03_filter_seurat.R \
--filter ${FILTER_NM} \
--seurat ${DATA_OUT}'unfiltered/seurat.rds'
--output ${FILTER_DATA_OUT} \

# Run integration
Rscript ${SCRIPTS}04_run_integration.R \
--input ${FILTER_DATA_OUT}filtered_seurat.rds \
--output ${FILTER_DATA_OUT} \
--figures ${RESULTS} \
--integration ${INTEGRATION_METHOD}

# Cluster barplots
Rscript ${SCRIPTS}05_barplots.R \
--metadata ${FILTER_DATA_OUT}integrated_seurat.csv.gz \
--figures ${RESULTS}cluster_barplots 

# Clone tables
Rscript ${SCRIPTS}05_tables.R \
--metadata ${FILTER_DATA_OUT}integrated_seurat.csv.gz \
--figures ${RESULTS}clone_tables \
--clones = ${WORKDIR}scripts/clones_of_interest.R"))




# 
# # Custom UMAPs
# Rscript ${SCRIPTS}08_custom_umaps.R \
# --metadata ${FILTER_DATA_OUT}integrated_seurat.csv.gz \
# --figures ${RESULTS}umap
# 
# # Differential expression
# Rscript ${SCRIPTS}09_differential_expression.R \
# --input ${FILTER_DATA_OUT}integrated_seurat.rds \
# --output ${RESULTS}differential_expression
