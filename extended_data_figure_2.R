#R code for extended data figure 2
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

#load scRNAseq seurat object, result of scRNAseq_analysis.R
mesentery <- readRDS("path_to_output_folder/CD45+/cd45+_cd_mesentery.rds")
#load the cluster markers, a result from scRNAseq_analysis.R
cell_markers <- read.csv("path_to_output_folder/CD45+/CD45+_CD_mesentery_marker_genes.csv")
#load tcr output csv file, result of vdj_clonality.py
tcr <- read.csv("path_to_folder/tcr_output.csv")

#sfig_2a
mesentery$sample_2 <- factor(mesentery$sample, 
                             levels = c("CD003_NI_MES", "CD003_INF_MES", "CD004_NI_MES", "CD004_INF_MES", "CD005_NI_MES", "CD005_INF_MES"),
                             labels = c("CD_03\nNon-inflamed", "CD_03\nInflamed", 
                                        "CD_04\nNon-inflamed", "CD_04\nInflamed",
                                        "CD_05\nNon-inflamed", "CD_05\nInflamed"))
sfig_2a <- DimPlot(mesentery, 
                   pt.size = 0.2, 
                   group.by = "sample_2") +
           theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 plot.title = element_blank(),
                 legend.text = element_text(size = 7),
                 legend.key.size = unit(5, "points"),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

#sfig_2b
mesentery$state_2 <- factor(mesentery$state, 
                           levels = c("CD Non-Inflamed", "CD Inflamed"),
                           labels = c("CD Non-inflamed", "CD Inflamed"))
sfig_2b <- DimPlot(mesentery,
                   pt.size = 0.2,
                   group.by = "state_2",
                   cols = c("#0000FF", "#FB0207")) +
           theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 plot.title = element_blank(),
                 legend.text = element_text(size = 7),
                 legend.key.size = unit(5, "points"),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

#sfig_2c
top_5 <- group_by(cell_markers, cluster, .add = TRUE) %>% slice(1:5) %>% select(cluster, gene)
mes <- NormalizeData(object = mesentery)
mes <- FindVariableFeatures(object = mes, nfeatures = 5000)
mes <- ScaleData(object = mes)
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

Idents(mesentery) <- "detailed_names"
sfig_2c <- DoHeatmap(subset(mes, downsample = 100), 
                     features = top_5$gene,
                     label = FALSE,
                     raster = FALSE) +
           theme(axis.text.y = element_text(size = 7),
                 legend.key.size = unit(5, "pt"),
                 legend.title = element_text(size = 7),
                 legend.text = element_text(size = 7),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
#sfig_2
make_clonal_graphs <- function(type, vdj) {
    data <- filter(vdj, grepl(type, metadata)) %>% group_by(patient) %>% count(is_clonal)
    data$is_clonal <- factor(data$is_clonal, labels = c("Non-Clonal", "Clonal"))
    data$patient <- factor(data$patient, 
                           levels = c("CTL_01", "CTL_02", "CD_03", "CD_04", "CD_05"),
                           labels = c("CTL\n01", "CTL\n02", "CD\n03", "CD\n04", "CD\n05"))
    plot <- ggplot(data, 
                  aes(x = patient, 
                      y = n, 
                      fill = is_clonal)) +  
           geom_bar(position = "fill", stat = "identity", color = "black", linewidth = 0.2) + 
           scale_fill_discrete(breaks = c("Clonal", "Non-Clonal")) +
           ggtitle(stringr::str_to_sentence(gsub(" ", "\n", gsub("NI", "Non-inflamed", gsub("II", "Inflamed", x = type))))) + 
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
                 legend.position = "bottom",
                 panel.background = element_blank(),
                 panel.grid = element_blank())
}

all_types = c("NI ileum", "II ileum", "NI mesentery", "II mesentery")
tcr_all_graphs <- lapply(X = all_types, vdj = tcr, FUN = make_clonal_graphs)
tcr_layout <- "AB
               CD"
sfig_2d <- tcr_all_graphs[[1]] + tcr_all_graphs[[2]] +
           tcr_all_graphs[[3]] + tcr_all_graphs[[4]] +
           plot_layout(design = tcr_layout, guides = "collect") & 
           theme(legend.position = "bottom",
                 legend.margin = margin(t = -10))

overall_layout <- "ACC
                   BCC
                   DCC"
supplemental_fig_2 <- sfig_2a + sfig_2b + wrap_elements(sfig_2c) + wrap_elements(sfig_2d) + plot_layout(design = overall_layout)

ggsave(filename = "path_to_output/supplemental_fig_1_r.pdf",
       plot = supplemental_fig_2,
       device = "pdf",
       units = "in",
       width = 7.5,
       height = 7.5,
       dpi = "print")