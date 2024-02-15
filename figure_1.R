#scRNAseq data for figure 1
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

#load scRNAseq seurat object, result of scRNAseq_analysis.R
mesentery <- readRDS("path_to_output_folder/CD45+/cd45+_cd_mesentery.rds")
#load bcr output csv file, result of vdj_clonality.py
bcr <- read.csv("path_to_folder/bcr_output.csv")

#annotation of the clusters
cluster_names <- c("Memory CD8+ T cells 1", "S100A8/S100A9 macrophages", "Memory CD4+ T cells 1", "Naïve B cells 1", 
                   "Memory B cells", "Memory CD8+ T cells 2", "Memory CD4+ T cells 2", "cDC2", "T cells", "Plasma cells",
                   "Naïve B cells 2", "NK/NKT cells", "Proliferating cells", "CD16+ macrophages", "CX3CR1+ T cells",
                   "LYVE-1+ macrophages", "Migrating DCs", "pDCs", "cDC1s")
Idents(mesentery) <- "seurat_clusters"
mesentery <- RenameIdents(mesentery, "0" = "0: Memory CD8+ T cells 1",
                                     "1" = "1: S100A8/S100A9 macrophages",
                                     "2" = "2: Memory CD4+ T cells 1",
                                     "3" = "3: Naïve B cells 1",
                                     "4" = "4: Memory B cells",
                                     "5" = "5: Memory CD8+ T cells 2",
                                     "6" = "6: Memory CD4+ T cells 2",
                                     "7" = "7: cDC2",
                                     "8" = "8: T cells",
                                     "9" = "9: Plasma cells",
                                     "10" = "10: Naïve B cells 2",
                                     "11" = "11: NK/NKT cells",
                                     "12" = "12: Proliferating cells",
                                     "13" = "13: CD16+ macrophages",
                                     "14" = "14: CX3CR1+ T cells",
                                     "15" = "15: LYVE-1+ macrophages",
                                     "16" = "16: Migrating DCs",
                                     "17" = "17: pDCs",
                                     "18" = "18: cDC1s")
mesentery$detailed_names <- mesentery@active.ident

#get the number of samples per cluster, seperated by state (inflamed vs non-inflamed)
stacked_bar_data <- mesentery@meta.data %>% group_by(state) %>% count(detailed_names)

my_colors <- scales::hue_pal()(19)

Idents(mesentery) <- "seurat_clusters"
plot_1c <- DimPlot(mesentery, pt.size = 0.2, label = TRUE, label.size = 2, label.box = TRUE) + 
           NoLegend() +
           theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

plot_1d <- ggplot(stacked_bar_data, aes(fill = detailed_names, 
                                        x = factor(state, level = c("CD Non-Inflamed", "CD Inflamed")), 
                                        y = n)) +  
           geom_bar(position = "fill", stat = "identity", color = "black", linewidth = 0.3) + 
           scale_fill_manual(values = c(my_colors)) +
           scale_x_discrete(labels = c("Non-inflamed\nMesentery", "Inflamed\nMesentery")) +
           theme(axis.text.x = element_text(size = 7, vjust = 5, face = "bold"),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(),
                 legend.key.size = unit(7, "points"),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 7),
                 legend.box.spacing = unit(-10, "points"),
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

#1e
identity_gene_list <- c("CD79A", "JCHAIN", "IGHA2", "CD3E", "CD4", "CD8A", "LYZ", "CCR7", "CLEC9A")
identity_plots <- lapply(identity_gene_list, FUN = function(x){
    FeaturePlot(mesentery, features = x, pt.size = 0.01) + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 7, face = "italic", vjust = -4),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.text = element_text(size = 6),
          legend.key.size = unit(5, "pt"))
})
plot_1e <- wrap_plots(identity_plots)

#1f
plot_1f <- DotPlot(mesentery, 
              features = c("TNF", "LTA", "LTB", "TNFSF14", "TNFRSF1A", "TNFRSF1B", "LTBR", "TNFRSF14"), 
              scale.by = "size",
              dot.scale = 3) + 
              scale_y_discrete(labels = cluster_names) + 
              scale_x_discrete(labels = c("TNF", "LTa", "LTb", "LIGHT", "TNFR1", "TNFR2", "LTbR", "HVEM")) + 
              scale_color_distiller(palette = "Blues", direction = 0) +
              theme(axis.title = element_blank(),
                    axis.text = element_text(size = 7),
                    axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(r = 7), face = "italic"),
                    legend.key.size = unit(5, "pt"),
                    legend.title = element_text(size = 7),
                    legend.text = element_text(size = 7),
                    plot.margin = margin(t = 0, r = 0, b = 0, l = 0))


#1g
make_clonal_graphs <- function(type, vdj) {
    data <- filter(vdj, grepl(type, metadata)) %>% group_by(patient) %>% count(is_clonal)
    data$is_clonal <- factor(data$is_clonal, labels = c("Non-Clonal", "Clonal"))
    data$patient <- factor(data$patient, 
                           levels = c("CTL_04", "CTL_05", "CD_01", "CD_02", "CD_03"),
                           labels = c("CTL\n01", "CTL\n02", "CD\n03", "CD\n04", "CD\n05"))
    plot<- ggplot(data, 
                  aes(x = patient, 
                      y = n, 
                      fill = is_clonal)) +  
           geom_bar(position = "fill", stat = "identity", color = "black", linewidth = 0.2) + 
           scale_fill_discrete(breaks = c("Clonal", "Non-Clonal")) +
           ggtitle(gsub(" ", "\n", type)) + 
           theme(plot.title = element_text(hjust = 0.5, vjust = -2, face = "bold", size = 7),
                 axis.text.x = element_text(size = 7, vjust = 4),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(),
                 legend.title = element_blank(),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
                 legend.key.size = unit(5, "points"),
                 legend.text = element_text(size = 6),
                 legend.box.spacing = unit(-10, "points"),
                 legend.position = "none",
                 panel.background = element_blank(),
                 panel.grid = element_blank())
}

all_types = c("Non-inflamed Ileum", "Inflamed Ileum", "Non-inflamed Mesentery", "Inflamed Mesentery")
bcr_all_graphs <- lapply(X = all_types, vdj = bcr, FUN = make_clonal_graphs)

#1h
isotype_graphs <- lapply(all_types, FUN = function(type) {
    data <- filter(bcr, grepl(type, metadata)) %>% group_by(patient) %>% count(igh_constant) %>% filter(igh_constant != "")
    data$patient <- factor(data$patient, 
                           levels = c("CTL_04", "CTL_05", "CD_01", "CD_02", "CD_03"),
                           labels = c("CTL\n01", "CTL\n02", "CD\n03", "CD\n04", "CD\n05"))
    data$igh_constant <- factor(data$igh_constant, levels = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHM"))
    plot<- ggplot(data, 
                  aes(x = patient, 
                      y = n, 
                      fill = igh_constant)) +  
           geom_bar(position = "fill", stat = "identity", color = "black", linewidth = 0.2) + 
           scale_fill_manual(values = c("IGHA1" = "#80b1d3",
                                        "IGHA2" = "#8dd3c7",
                                        "IGHG1" = "#ffffb3",
                                        "IGHG2" = "#bebada",
                                        "IGHG3" = "#fdb462",
                                        "IGHG4" = "#b3de69",
                                        "IGHD" = "#fccde5",
                                        "IGHM" = "#fb8072")) +
           ggtitle(gsub(" ", "\n", type)) + 
           theme(plot.title = element_text(hjust = 0.5, vjust = -2, face = "bold", size = 7),
                 axis.text.x = element_text(size = 7, vjust = 4),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 6),
                 legend.key.size = unit(5, "points"),
                 legend.box.spacing = unit(-10, "points"),
                 legend.position = "none",
                 panel.background = element_blank(),
                 panel.grid = element_blank())
})

bcr_layout <- "ABCD
               EFGH"
plot_1g_h <- bcr_all_graphs[[1]] + bcr_all_graphs[[2]] + isotype_graphs[[1]] + isotype_graphs[[2]] + 
             bcr_all_graphs[[3]] + bcr_all_graphs[[4]] + isotype_graphs[[3]] + isotype_graphs[[4]] + 
             plot_layout(design = bcr_layout)

fig_1 <- (plot_spacer() / wrap_elements(full = plot_1c) / plot_1d) | (wrap_elements(full = plot_1e) / plot_1f / wrap_elements(full = plot_1g_h))

ggsave(filename = "path_to_output_folder/fig_1_r.pdf",
       plot = fig_1,
       device = "pdf",
       units = "in",
       width = 7.2,
       height = 7.5,
       dpi = "print")