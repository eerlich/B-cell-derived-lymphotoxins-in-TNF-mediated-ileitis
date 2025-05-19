library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)

setwd("path_to_output_dir/")

human <- readRDS("cd45+_cd_only_mesentery.rds")
mouse <- readRDS("mln_lt_blocking.rds")

human <- RenameIdents(human, "0" = "Memory CD8+ T cells 1",
                             "1" = "S100A8/S100A9 macrophages",
                             "2" = "Memory CD4+ T cells 1",
                             "3" = "Naïve B cells 1",
                             "4" = "Memory B cells",
                             "5" = "Memory CD8+ T cells 2",
                             "6" = "Memory CD4+ T cells 2",
                             "7" = "cDC2",
                             "8" = "T cells",
                             "9" = "Plasma cells",
                             "10" = "Naïve B cells 2",
                             "11" = "NK/NKT cells",
                             "12" = "Proliferating cells",
                             "13" = "CD16+ macrophages",
                             "14" = "CX3CR1+ T cells",
                             "15" = "LYVE-1+ macrophages",
                             "16" = "Migrating DCs",
                             "17" = "pDCs",
                             "18" = "cDC1s")
human$detailed_names <- human@active.ident

#supplemental figure 5b
sfig5_b <- DimPlot(human,
                   reduction = "umap",
                   group.by = "state",
                   cols = c("#ef8a62", "#67a9cf")) +
           coord_fixed() +
           theme(plot.title = element_blank(),
                 legend.text = element_text(size = 5),
                 legend.position= c(.7, .8),
                 axis.line = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text = element_blank(),
                 legend.key.size = unit(6, "pt"),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

#supplemental figure 5c
Idents(human) <- "seurat_clusters"
sfig5_c <- DimPlot(human,
                   reduction = "umap",
                   pt.size = 0.2, 
                   label = TRUE, 
                   label.size = 2, 
                   label.box = TRUE) +
           coord_fixed() +
           scale_color_discrete(labels = cluster_names_num) +
           guides(color = guide_legend(override.aes = list(size = 2))) +
           theme(plot.title = element_blank(),
                 legend.text = element_text(size = 5),
                 axis.line = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text = element_blank(),
                 legend.key.size = unit(0.1, "cm"),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

#supplemental figure 5d
human_markers <- read.csv("/Users/eerlich/Desktop/projects/phd/CD45+_CD_only_mesentery_marker_genes.csv")
top_3 <- group_by(human_markers, cluster, .add = TRUE) |>
         slice(1:3) |> 
         select(cluster, gene)
Idents(human) <- "detailed_names"
sfig5_d <- DoHeatmap(subset(human, downsample = 100), 
                     features = top_3$gene,
                     label = FALSE,
                     raster = FALSE) +
           scale_fill_viridis() +
           theme(axis.text.y = element_text(size = 5),
                 legend.key.size = unit(5, "pt"),
                 legend.title = element_text(size = 7),
                 legend.text = element_text(size = 7),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

#supplemental figure 5e
top_genes <- c("Cd79a", "Cd79b", "Fcmr", "Cd19", "Ighm",
               "Lat", "Cd3e", "Tcf7", "Il7r", "Cd4",
               "Arhgap26", "Wdfy4", "Pax5", "Myo1e", "Chst3",
               "Cd8b1", "Cd8a", "Lef1", "Nkg7", "Thy1",
               "Themis", "Bcl11b", "Itk", "Cd247", "Tox",
               "Fscn1", "Basp1", "Cxcl16", "Anxa3", "Cd63",
               "Ms4a4b", "Ccl5", "Il2rb", "Ctsw", "Cxcr3", 
               "Fcer1g", "Wfdc17", "Mpeg1", "Lyz2", "Csf1r")
mouse <- ScaleData(mouse, features = top_genes)
sfig5_e <- DoHeatmap(subset(mouse, downsample = 100), 
                     features = top_genes,
                     label = FALSE,
                     raster = FALSE) +
           scale_fill_viridis() +
           theme(axis.text.y = element_text(size = 5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
                 legend.key.size = unit(5, "pt"),
                 legend.title = element_text(size = 7),
                 legend.text = element_text(size = 7),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

#supplemental figure 5f
sfig5_f <- FeatureScatter(subset(mouse, ident = "B cells 1"), 
                          feature1 = "Lta", 
                          feature2 = "Tnf") +
           NoLegend() +
           coord_fixed() +
           ggtitle("B cell cluster 1") +
           theme(plot.title = element_text(size = 8),
                 axis.text = element_text(size = 5),
                 axis.title = element_text(size = 7, face = "italic"),
                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

ggsave(filename = "sfig5_b_r.pdf",
       plot = sfig5_b,
       device = "pdf",
       units = "in",
       width = 2,
       height = 2,
       dpi = "print")

ggsave(filename = "sfig5_c_r.pdf",
       plot = sfig5_c,
       device = "pdf",
       units = "in",
       width = 2,
       height = 2,
       dpi = "print")

ggsave(filename = "sfig5_d_r.pdf",
       plot = sfig5_d,
       device = "pdf",
       units = "in",
       width = 3,
       height = 4,
       dpi = "print")

ggsave(filename = "sfig5_e_r.pdf",
       plot = sfige_f,
       device = "pdf",
       units = "in",
       width = 2.75,
       height = 2.75,
       dpi = "print")

ggsave(filename = "sfig5_f_r.pdf",
       plot = sfig5_f,
       device = "pdf",
       units = "in",
       width = 2,
       height = 2,
       dpi = "print")