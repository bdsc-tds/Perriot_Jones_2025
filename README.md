# Perriot and Jones, 2025
### Neuron-reactive KIR+ CD8+ T cells display an encephalitogenic transcriptional program in autoimmune encephalitis

Code for R analyses for Perriot and Jones et al, 2025.  The TCR-seq data used in the study is available at NCBI GEO with accession [GSE263666](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263666).  Running the analyses involving tested clones requires a list of the TCR sequences called "merged_clone_sets.csv" and a script containing a clone of interest "clones_of_interest.R" which we do not include in this github repository.  Please contact the corresponding authors for access to this table.  

The pre-processing scripts prefixed with numbers must be run in order.  The 
other scripts are stand-alone (assuming the pre-processing has been completed).

The starting point for the R analyses are the cellbender adjusted counts tables.
Cellbender was installed in an anaconda environment, and R analyses were run in
a singularity container.  The yaml file to regenerate the conda environment and
the container definition are found in the "envs" folder.  
