#khroma::scale_color_discreterainbow(reverse = TRUE)

#p1_rem <- DimPlot(ri01, group.by = "Sample",
#                  cols = c("lightgray", "#DC050C"), order = "Ri01_5m") +
#    guides(color = "none")


# Minimal feature plot, yellow red theme
#FeaturePlot(Ri01_dis, features = fig2b_markers, keep.scale = "all",
#            cols = c("lightgray", "gold", "firebrick")) +
#    plot_layout(guides = "collect") &
#    theme(axis.line = element_blank(),
#          axis.title = element_blank(),
#          axis.text = element_blank(),
#          axis.ticks = element_blank()) &
#    labs(color = "log2 Exp")


# Get data 
#av_expr <- AggregateExpression(Ri01_dis, features = fig4c_markers,
#                               return.seurat = TRUE)
#av_expr <- FetchData(av_expr, vars = fig4c_markers)
#av_expr <- t(av_expr)