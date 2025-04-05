# Perriot_Jones_2025
Code for R analyses for Perriot and Jones et al, 2025.  The TCR-seq data used in the study will be made availabe upon acceptance of the manuscript.  Running the analyses involving tested clones requires a list of the TCR sequences which we do not include in this github repository.  Please contact the corresponding authors for access to this table.  

Currently, this repository includes exploratory code.  We are in the process of 
tidying the repository to just include the analyses in the submitted manuscript,
and renaming scripts to match the figures in the manuscript. 

The pre-processing scripts prefixed with numbers must be run in order.  The 
other scripts are stand-alone (assuming the pre-processing has been completed).

The starting point for the R analyses are the cellbender adjusted counts tables.
Cellbender was installed in an anaconda environment, and R analyses were run in
a singularity container.  The yaml file to regenerate the conda environment and
the container definition are found in the "envs" folder.  
