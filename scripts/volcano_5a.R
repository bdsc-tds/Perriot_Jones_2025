# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ComplexHeatmap")
library("ggplot2")
library("tidyverse")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--de',  
                    help = 'csv with differential expression results')
parser$add_argument('--workdir',
                    help = "Working directory, for loading scripts")

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/funcs_EnhancedVolcano_mod.R"))

# ----------------------------------------------------------------------------
# main ----
markers <- function(){
    return(c("IKZF2", "KLRC2", "KLRC3", "KLRK1", "KIR2DL1", "KIR2DL3", "IL7R",
             "ZAP70", "FOS", "JUN", "RELB", "EGR1", "EGR2",
             "EGR3", "JUND", "FOSB", "INPP5D", "CBLB",
             "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5",
             "CD69", "CD74", "CD44", "VCAM1", "ICAM1",
             "TNF", "IFNG", "GZMA", "GZMH", "CCL4", "CCL3", "CSF1",
             "ZFP36L1", "NFKBIA"))
}

main <- function(args){
    de <- read_csv(file.path("results/one_tcr_beta/clone_analyses/",
                             "reactive_Ri_v_HD/diff_expr_full.csv"))

    pdf("results/one_tcr_beta/figures_in_paper/volcano_fig_5/volcano_reactive_Ri_v_HD.pdf")
    p <- volc_mod(de,
                  labels = markers(),
                  labelsFrom = "gene",
                  labSize = 5.5,
                  force = 50,
                  col = c("grey50", "red2","royalblue")) +
        theme(axis.title = element_text(size = 30),
              axis.text = element_text(size = 24))
    print(p)
    dev.off()
}

# ----------------------------------------------------------------------------
main(args)

# Fig 5a = reactive Ri v HD