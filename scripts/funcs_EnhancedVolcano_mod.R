# Modified from EnhancedVolcano for better control over labels
# note - drawConnectors determines whether geom_text_repel or geom_text is used
# Argument force in geom_text_repel controls overlaps,
# max.overlaps = Inf forces all requested labels to be shown
# min.segment.length = 0 forces all lines to be drawn
# default text is bold, use "plain" for un-bolded 
# TO DO - PROPERLY CONSTRUCT LABEL COLOUR
# ... args for ggrepel

library("ggrepel")

volc_mod <- function(tt, xlab, ylab, x = "logFC", y = "adj_pval",
                     FCcutoff = 0.25, pCutoff = 0.05, pointSize = 2, labSize = 5,
                     shape = 19, colAlpha = 0.5, legendIconSize = 5,
                     col = c("grey30", "red2","royalblue"),
                     #col = c("red2","royalblue"), # "#008000"
                     ...){
    
    # Set "Sig" for colours
    tt$Sig <- "ns"
    tt$Sig[(tt[[y]] < pCutoff) & tt[[x]] > 0] <- "up"
    tt$Sig[(tt[[y]] < pCutoff) & tt[[x]] < 0] <- "down"
    tt$Sig <- factor(tt$Sig, levels = c("ns", "up", "down"))
    
    pCutoff = max(tt[which(tt[y] <= pCutoff), y])
    
    p <- ggplot(tt, aes(x = !!sym(x), y = -log10(!!sym(y)))) + 
        volcano_theme() +
        guides(colour = guide_legend(override.aes = list(size = legendIconSize))) + 
        geom_point(aes(color = Sig), 
                   alpha = colAlpha, 
                   shape = shape,
                   size = pointSize,
                   na.rm = TRUE) + 
        scale_color_manual(values = col) +
        guides(color = "none") +
        geom_vline(xintercept = c(-FCcutoff, FCcutoff),
                   linetype = "dashed", colour = "darkgray", linewidth = 0.4) +
        geom_hline(yintercept = -log10(pCutoff), 
                   linetype = "dashed", colour = "darkgray", linewidth = 0.4) +
        labs(x = xlab, y = ylab)
    p
}
    

volcano_theme <- function(){
    th <- theme_bw(base_size = 24) +
        theme(axis.title = element_text(size = 18),
              legend.position = "top", 
              legend.key.size = unit(0.5, "cm"),
              legend.text = element_text(size = 14), 
              legend.title = element_blank(),
              panel.grid = element_blank(),
              axis.line = element_line(linewidth = 0.8, 
                                       colour = "black"),
              panel.border = element_blank(),
              panel.background = element_blank())
    return(th)
}
    