#CD45+ on human mesentery

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd("path_to_output_folder/CD45+/")

CD003_NI_MES_path <- "path_to/CD003_NI_MES/outs/per_sample_outs/CD003_NI_MES/count/sample_filtered_feature_bc_matrix/"
CD003_INF_MES_path <- "path_to/CD003_INF_MES/outs/per_sample_outs/CD003_INF_MES/count/sample_filtered_feature_bc_matrix/"
CD004_NI_MES_path <- "path_to/CD004_NI_MES/outs/per_sample_outs/CD004_NI_MES/count/sample_filtered_feature_bc_matrix/"
CD004_INF_MES_path <- "path_to/CD004_INF_MES/outs/per_sample_outs/CD004_INF_MES/count/sample_filtered_feature_bc_matrix/"
CD005_NI_MES_path <- "path_to/CD005_NI_MES/outs/per_sample_outs/CD005_NI_MES/count/sample_filtered_feature_bc_matrix/"
CD005_INF_MES_path <- "path_to/CD005_INF_MES/outs/per_sample_outs/CD005_INF_MES/count/sample_filtered_feature_bc_matrix/"

#list of file paths
file_paths <- c(CD003_NI_MES_path, CD003_INF_MES_path, CD004_NI_MES_path, CD004_INF_MES_path, CD005_NI_MES_path, CD005_INF_MES_path)

sample_list <- lapply(file_paths, FUN = function(file_path){
    sample_data <- Read10X(file_path)
    obj <- CreateSeuratObject(counts = sample_data, project = "Human Mesentery",
          min.cells = 3, min.features = 200)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    return(obj)
})

#annotate the samples
sample_list[[1]]$sample <- "CD003_NI_MES"
sample_list[[1]]$patient <- "CD003"
sample_list[[1]]$state <- "CD Non-Inflamed"
sample_list[[1]]$tissue <- "Mesentery"

sample_list[[2]]$sample <- "CD003_INF_MES"
sample_list[[2]]$patient <- "CD003"
sample_list[[2]]$state <- "CD Inflamed"
sample_list[[2]]$tissue <- "Mesentery"

sample_list[[3]]$sample <- "CD004_NI_MES"
sample_list[[3]]$patient <- "CD004"
sample_list[[3]]$state <- "CD Non-Inflamed"
sample_list[[3]]$tissue <- "Mesentery"

sample_list[[4]]$sample <- "CD004_INF_MES"
sample_list[[4]]$patient <- "CD004"
sample_list[[4]]$state <- "CD Inflamed"
sample_list[[4]]$tissue <- "Mesentery"

sample_list[[5]]$sample <- "CD005_NI_MES"
sample_list[[5]]$patient <- "CD005"
sample_list[[5]]$state <- "CD Non-Inflamed"
sample_list[[5]]$tissue <- "Mesentery"

sample_list[[6]]$sample <- "CD005_INF_MES"
sample_list[[6]]$patient <- "CD005"
sample_list[[6]]$state <- "CD Inflamed"
sample_list[[6]]$tissue <- "Mesentery"

#make QC plots
QC_plots <- lapply(sample_list, FUN = function(sample_obj){
    p1 <- VlnPlot(sample_obj, features = "nFeature_RNA") +
          geom_hline(yintercept= c(500, 5000), linetype = "dashed", color = "red") +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                legend.position = "none"
          )
    p2 <- VlnPlot(sample_obj, features = "nCount_RNA") +
          geom_hline(yintercept = c(25000), linetype = "dashed", color = "red") +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                legend.position = "none"
          )
    p3 <- VlnPlot(sample_obj, features = "percent.mt") +
          geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                legend.position = "none"
          )
    p4 <- ggplot(sample_obj@meta.data, aes(x=nCount_RNA)) + geom_histogram(binwidth = 100)
    p5 <- ggplot(sample_obj@meta.data, aes(x=nFeature_RNA)) + geom_histogram(binwidth = 100)
    total <- (p1 | p2 | p3) / (p4 | p5) + plot_annotation(title = sample_obj@meta.data$sample[1])
})

pdf("all_mesentery_quality_control_plots.pdf")
QC_plots
dev.off()

#perform QC, normalize, and identify variable features for every sample seperately
sample_list <- lapply(sample_list, FUN = function(x){
    x <- subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA < 25000 & percent.mt < 10)
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})

#perform integration
variable_features <- SelectIntegrationFeatures(object.list = sample_list)
cell_anchors <- FindIntegrationAnchors(object.list = sample_list, anchor.features = variable_features)
cells <- IntegrateData(anchorset = cell_anchors)
DefaultAssay(cells) <- "integrated"

#standard workflow for visualization and clustering
cells <- ScaleData(cells)
cells <- RunPCA(cells, npcs = 50, verbose = FALSE)

#plots for clustering
p1 <- ElbowPlot(cells, ndims = 50)
cells <- JackStraw(cells, num.replicate = 100, dims = 50)
cells <- ScoreJackStraw(cells, dims = 1:50)
p2 <- JackStrawPlot(cells, dims = 1:50)
pdf("mes_dimensionality_plots.pdf")
p1
p2
DimHeatmap(cells, dims = 1:12, cells = 500, balanced = TRUE)
DimHeatmap(cells, dims = 13:24, cells = 500, balanced = TRUE)
DimHeatmap(cells, dims = 25:36, cells = 500, balanced = TRUE)
DimHeatmap(cells, dims = 37:50, cells = 500, balanced = TRUE)
dev.off()

#cluster the data
cells <- RunUMAP(cells, reduction = "pca", dims = 1:33)
cells <- FindNeighbors(cells, reduction = "pca", dims = 1:33)
cells <- FindClusters(cells, resolution = 0.5)

#save the data
DefaultAssay(cells) <- "RNA"
saveRDS(cells, file = "cd_unfiltered_mesentery.rds")
cells_markers <- FindAllMarkers(cells, only.pos = T)
write.csv(cells_markers, "cd_unfiltered_mesentery_markers.csv")

#Visualize CD45- cells to crosscheck the results from the cell markers
pdf("cd_mesentery_cd45+_purity.pdf")
DimPlot(cells, label = T)
FeaturePlot(cells, features = c("PTPRC", "CD79A"), label = TRUE)
dev.off()

#remove the CD45- clusters
cells <- subset(cells, idents = c("3", "9", "15", "16"), invert = TRUE)

#resplit the samples
sample_list <- lapply(unique(cells@meta.data[["sample"]]), FUN = function(x) {
    y <- subset(cells, subset = sample == x)
    y <- NormalizeData(y, verbose = FALSE)
    y <- FindVariableFeatures(y, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    return(y)
})

#rejoin the samples
variable_features <- SelectIntegrationFeatures(object.list = sample_list)
cell_anchors <- FindIntegrationAnchors(object.list = sample_list, anchor.features = variable_features)
cells <- IntegrateData(anchorset = cell_anchors)
DefaultAssay(cells) <- "integrated"

#redo the workflow for visualization and clustering
cells <- ScaleData(cells, verbose = FALSE)
cells <- RunPCA(cells, npcs = 50, verbose = FALSE)

#plots for clustering
p1 <- ElbowPlot(cells, ndims = 50)
cells <- JackStraw(cells, num.replicate = 100, dims = 50)
cells <- ScoreJackStraw(cells, dims = 1:50)
p2 <- JackStrawPlot(cells, dims = 1:50)
pdf("CD45+_cd_mesentery_Dimensionality_plots.pdf")
p1
p2
DimHeatmap(cells, dims = 1:12, cells = 500, balanced = TRUE)
DimHeatmap(cells, dims = 13:24, cells = 500, balanced = TRUE)
DimHeatmap(cells, dims = 25:36, cells = 500, balanced = TRUE)
DimHeatmap(cells, dims = 37:50, cells = 500, balanced = TRUE)
dev.off()

#cluster the data
cells <- RunUMAP(cells, reduction = "pca", dims = 1:19)
cells <- FindNeighbors(cells, reduction = "pca", dims = 1:19)

pdf("CD45+_CD_Mesentery_Cluster_Resolution_plots.pdf")
lapply(seq(from = 0.1, to = 1.0, by = 0.1), FUN = function(res){
      DefaultAssay(cells) <- "integrated"
      cells <- FindClusters(cells, resolution = res)
      p1 <- DimPlot(cells, pt.size = 1.5, label = TRUE, label.size = 10) + labs(title = paste("Cluster Resolution:", res))
      return(p1)
})
dev.off()

#get cluster markers
DefaultAssay(cells) <- "integrated"
cells <- FindClusters(cells, resolution = 0.6)
DefaultAssay(cells) <- "RNA"
cell_markers <- FindAllMarkers(cells, only.pos = TRUE)
cell_markers <- filter(cell_markers, p_val_adj < 0.05)
write.csv(cell_markers, file = "CD45+_CD_mesentery_marker_genes.csv", row.names = FALSE)

saveRDS(cells, file = "cd45+_cd_mesentery.rds")