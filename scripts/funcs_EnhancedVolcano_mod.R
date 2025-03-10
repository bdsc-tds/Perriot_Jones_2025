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
                     FCcutoff = 1, pCutoff = 0.05, pointSize = 2, labSize = 5,
                     shape = 19, colAlpha = 0.5, legendIconSize = 5,
                     col = c("red2","royalblue"), # "#008000"
                     labels = NULL, labelsFrom = "name", labelColour = NULL,
                     force = 5, fontface = "bold", ...){
    
    # Set "Sig" for colours
    tt$Sig <- "ns"
    tt$Sig[(tt[[y]] < pCutoff) & tt[[x]] > 0] <- "up"
    tt$Sig[(tt[[y]] < pCutoff) & tt[[x]] < 0] <- "down"
    tt$Sig <- factor(tt$Sig, levels = c("ns", "up", "down"))
    
    pCutoff = max(tt[which(tt[y] <= pCutoff), y])
    
    p <- ggplot(tt, aes(x = !!sym(x), y = -log10(!!sym(y)))) + 
        volcano_theme() +
        guides(colour = guide_legend(override.aes = list(size = legendIconSize))) + 
        geom_point(color = "grey40", 
                   alpha = colAlpha, 
                   shape = shape,
                   size = pointSize,
                   na.rm = TRUE) + 
        geom_vline(xintercept = c(-FCcutoff, FCcutoff),
                   linetype = "dashed", colour = "darkgray", linewidth = 0.4) +
        geom_hline(yintercept = -log10(pCutoff), 
                   linetype = "dashed", colour = "darkgray", linewidth = 0.4) +
        labs(x = xlab, y = ylab)
    
    if (! is.null(labels)){
        
        tt_subs <- tt[tt[[labelsFrom]] %in% labels, ]
        
        # Add coloured points ----
        p <- p + geom_point(data = tt_subs, aes(color = Sig),
                            size = pointSize + 0.5) +
            scale_color_manual(values = col) + #labels = legendLabels)  +
            guides(color = "none") 
        
        # Add repel labels ----    
        tt_up <- tt[tt[[labelsFrom]] %in% labels & tt$Sig == "up", ]
        tt_down <- tt[tt[[labelsFrom]] %in% labels & tt$Sig == "down", ]
        
        # Do not allow labels to start underneath p-val cutoff
        x_range <- range(tt[[x]])
        y_range <- c(-log10(pCutoff), max(na.omit(-log10(tt[[y]]))))

        repel_args <- list(size = labSize,
                           aes(label = !! sym(labelsFrom)),
                           force = force,
                           max.overlaps = Inf,
                           min.segment.length = 0,
                           fontface = fontface)
        p <- p +
            do.call(ggrepel::geom_text_repel,
                    utils::modifyList(repel_args,
                                      c(list(data = tt_up,
                                             xlim = c(FCcutoff, x_range[2]),
                                             ylim = y_range),
                                        list(...)))) +
            do.call(ggrepel::geom_text_repel,
                    utils::modifyList(repel_args,
                                      c(list(data = tt_down,
                                             xlim = c(x_range[1], -FCcutoff),
                                             ylim = y_range),
                                        list(...))))
    }
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
    