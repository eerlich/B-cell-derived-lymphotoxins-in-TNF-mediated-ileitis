library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd("path_to_output_dir/")

human <- readRDS("cd45+_cd_only_mesentery.rds")
mouse <- readRDS("mln_lt_blocking.rds")

#annotation of the human clusters
cluster_names <- c("Memory CD8+ T cells 1", "S100A8/S100A9 macrophages", "Memory CD4+ T cells 1", "Naïve B cells 1", 
                   "Memory B cells", "Memory CD8+ T cells 2", "Memory CD4+ T cells 2", "cDC2", "T cells", "Plasma cells",
                   "Naïve B cells 2", "NK/NKT cells", "Proliferating cells", "CD16+ macrophages", "CX3CR1+ T cells",
                   "LYVE-1+ macrophages", "Migrating DCs", "pDCs", "cDC1s")

#figure 7a
fig7_a <- DotPlot(human, 
                  features = c("TNF", "LTA", "LTB", "TNFSF14", "TNFRSF1A", "TNFRSF1B", "LTBR", "TNFRSF14"), 
                  scale.by = "size",
                  cols = "RdBu",
                  dot.scale = 3) + 
          scale_y_discrete(labels = cluster_names) + 
          scale_x_discrete(labels = c("TNF", "LT\u03B1", "LT\u03B2", "LIGHT", "TNFR1", "TNFR2", "LT\u03B2R", "HVEM")) +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 7),
                axis.text.x = element_text(angle = 45, hjust = 1.1, margin = margin(r = 7), face = "italic"),
                legend.key.size = unit(5, "pt"),
                legend.title = element_text(size = 7),
                legend.text = element_text(size = 7),
                plot.margin = margin(t = 0, r = 0, b = 0, l = 0))


#figure 7c
fig7_c_gene_list <- c("Ighd", "Cd3e", "Fscn1", "Ltb", "Lta", "Tnf")

fig7_c_list <- lapply(fig7_c_gene_list, FUN = function(x){
    p <- FeaturePlot(mouse,
                    features = x, 
                    reduction = "umap", 
                    cols = c("#ebebeb", "#cb181d")) +
         guides(col = guide_colourbar(barheight = unit(0.4, "in"), 
                                      barwidth = unit(0.05, "in"))) +
         theme(plot.title = element_text(face = "italic", size = 10),
               legend.text = element_text(size = 7),
               axis.line = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank(),
               axis.text = element_blank(),
               plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
    p$layers[[1]]$aes_params$shape <- 16
    p$layers[[1]]$aes_params$size <- 0.001
    return(p)
})

fig7_c <- wrap_plots(fig7_c_list, ncol = 2)

#figure 7d
fig7_d <- DimPlot(mouse,
                  reduction = "umap",
                  split.by = "sample") &
          coord_fixed() &
          theme(strip.text.x = element_text(size = 6),
               legend.text = element_text(size = 5),
               axis.line = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank(),
               axis.text = element_blank(),
               legend.key.size = unit(6, "pt"),
               plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

#figure 7e
#subset the clusters that are high in LTa expression
subset_cells <- subset(mouse, idents = c("B cells 1", "B cells 2", "CD4 T cells 1"))

subset_cells$cluster_names <- factor(subset_cells$cluster_names, 
                                     levels = c("B cells 1", "B cells 2", "CD4 T cells 1"))

Idents(subset_cells) <- "cluster_names"

#y_axis labels
y_axis_names <- c("Tnf\u0394ARE Anti-LT\u03B13", "Tnf\u0394ARE Isotype", "WT Anti-LT\u03B13", "WT Isotype",
                  "Tnf\u0394ARE Anti-LT\u03B13", "Tnf\u0394ARE Isotype", "WT Anti-LT\u03B13", "WT Isotype",
                  "Tnf\u0394ARE Anti-LT\u03B13", "Tnf\u0394ARE Isotype", "WT Anti-LT\u03B13", "WT Isotype")

fig7_e <- DotPlot(subset_cells, 
                  features = c("Ltb", "Lta", "Tnf", "Tnfrsf1a", "Tnfrsf1b", "Ltbr", "Tnfrsf14",
                               "Igha", "Il4i1", "Fcer2a", "Cd83", "Bhlhe40", "C1qbp", "Nr4a1"), 
                  scale.by = "size",
                  split.by = "sample",
                  cols = "RdBu",
                  dot.scale = 3) +
          scale_y_discrete(labels = y_axis_names) + 
          scale_x_discrete(labels = c("Lt\u03B2", "Lt\u03B1", "Tnf", "Tnfr1", "Tnfr2", 
                                      "Lt\u03B2r", "Light", "Igha", "Il4i1", "Fcer2a", 
                                      "Cd83", "Bhlhe40", "C1qbp", "Nr4a1")) +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 7),
                axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(r = 7), face = "italic"),
                legend.key.size = unit(5, "pt"),
                legend.title = element_text(size = 7),
                legend.text = element_text(size = 7),
                plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

ggsave(filename = "fig_7c_r.pdf",
       plot = fig7_c,
       device = cairo_pdf,
       units = "in",
       width = 2.6,
       height = 2.4,
       dpi = "print")

ggsave(filename = "fig_7_d_r.pdf",
       plot = fig7_d,
       device = cairo_pdf,
       units = "in",
       width = 4.4,
       height = 2,
       dpi = "print")

ggsave(filename = "fig_7_e_r.pdf",
       plot = fig7_e,
       device = cairo_pdf,
       units = "in",
       width = 4.4,
       height = 2,
       dpi = "print")