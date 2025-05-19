#scRNAseq analysis on just the MLN samples

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(reticulate)
library(sceasy)

setwd("path_to_output_dir/MLN")

#load in the 10X sample_filtered_feature_bc_matrix
input_data_dir <- file.path("..", "processed_data")
sample_name_dir_list <- list.files(input_data_dir)

#loop through the different 10X datasets
sample_list <- c()
sample_name_list <- c()
metadata_list <- c()
count <- 1
for (i in 1:length(sample_name_dir_list)) {
    sample_name <- str_split(sample_name_dir_list[[i]], "_", n = 2)[[1]][2]
    #set up the metadata for each sample
    temp <- str_split_1(sample_name, "_")
    if (temp[[2]] == "LN") {
        data <- Read10X(file.path(input_data_dir, sample_name_dir_list[[i]], "counts", "sample_filtered_feature_bc_matrix"))

        #make cell names unique for every sample
        data@Dimnames[[2]] <- paste0(data@Dimnames[[2]], "_", count)

        meta <- data.frame("cell_names" = data@Dimnames[[2]])
        meta <- meta %>%
                mutate(sample = sample_name,
                       donor = temp[[1]],
                       tissue = temp[[2]],
                       treatment = temp[[3]])

        sample_list[[count]] <- data
        sample_name_list[[count]] <- sample_name
        metadata_list[[count]] <- meta
        count = count + 1
    }
}

metadata <- bind_rows(metadata_list)
metadata <- tibble::column_to_rownames(metadata, var = "cell_names")
names(sample_list) <- sample_name_list

#make one object with the right metadata
data <- CreateSeuratObject(counts = sample_list, 
                           meta.data = metadata, 
                           project = "LTa blocking", 
                           min.cells = 3, 
                           min.features = 200)

#add the percentage of mitochrondria reads
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")

#run scDblFinder
data <- NormalizeData(data, verbose = FALSE)
data <- FindVariableFeatures(data, verbose = FALSE)
data <- ScaleData(data, verbose = FALSE)
data <- RunPCA(data, verbose = FALSE)
data <- IntegrateLayers(object = data,
                        method = CCAIntegration,
                        orig.reduction = "pca",
                        new.reduction = "integrated_cca")
data <- FindNeighbors(data, dims = 1:30, reduction = "integrated_cca")
data <- RunUMAP(data, dims = 1:30, reduction = "pca", reduction.name = "unfiltered_mln_umap")
data <- FindClusters(data, resolution = 1, cluster.name = "unfiltered_mln_clusters")
saveRDS(data, "mln_unfiltered.rds")

#convert the Seurat object to SingleCellExperiment using sceasy, save to file, then reload
data <- JoinLayers(data)
convertFormat(data, from = "seurat", to = "sce", outFile = "mln_unfiltered_sce.rds")
rm(list = ls())
gc()
sce <- readRDS("mln_unfiltered_sce.rds")

#doublet detection using scDblFinder
library(scDblFinder)
#in case the doublet estimate was off, want it to be called experimentally
sce <- scDblFinder(sce, samples = "sample", clusters = "unfiltered_mln_clusters", dbr.sd = 1, verbose = TRUE)
saveRDS(sce, "mln_unfiltered_sce.rds")

#add the doublet metadata to the seurat object
data <- readRDS("mln_unfiltered.rds")
data$sc_dbl_score <- sce$scDblFinder.score
data$sc_dbl_class <- sce$scDblFinder.class
saveRDS(data, "mln_unfiltered.rds")

#restart R
quit()
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(RColorBrewer)

#set working directory
setwd("path_to_output_dir/MLN")
data <- readRDS("mln_unfiltered.rds")

#add in cell cycle scoring
s_genes <- str_to_sentence(cc.genes$s.genes)
g2m_genes <- str_to_sentence(cc.genes$g2m.genes)
data <- JoinLayers(data)
data <- CellCycleScoring(data, s.features = s_genes, g2m.features = g2m_genes)

saveRDS(data, "mln_unfiltered.rds")


#filter out the low quality cells
data <- subset(data,
               subset = nFeature_RNA > 750 &
               nFeature_RNA < 4500 &
               nCount_RNA < 20000 &
               percent.mt < 5)

#remove clusters >50% doublets
data <- subset(data,
               ident = c("19", "24"),
               invert = TRUE)

#remove doublets individually
data <- subset(data,
               subset = sc_dbl_class == "singlet")

saveRDS(data, "mln_lt_blocking.rds")

#recluster without the suboptimal cells
data <- NormalizeData(data, verbose = FALSE)
data <- FindVariableFeatures(data, verbose = FALSE)
data <- ScaleData(data, verbose = FALSE)
data <- RunPCA(data, verbose = FALSE)
data <- IntegrateLayers(object = data,
                        method = CCAIntegration,
                        orig.reduction = "pca",
                        new.reduction = "integrated_cca")
ElbowPlot(data, ndims = 30)
data <- FindNeighbors(data, dims = 1:25, reduction = "integrated_cca")
data <- RunUMAP(data, dims = 1:25, reduction = "pca", reduction.name = "umap")
data <- FindClusters(data, resolution = 0.1, cluster.name = "mln_clusters")
data <- JoinLayers(data)

data <- RenameIdents(data, "0" = "B cells 1",
                           "1" = "CD4 T cells 1",
                           "2" = "B cells 2",
                           "3" = "CD8 T cells",
                           "4" = "CD4 T cells 2",
                           "5" = "Dendritic cells",
                           "6" = "Activated T cells",
                           "7" = "Monocytes")
data[["cluster_names"]] <- Idents(data)

data$treatment <- factor(data$treatment, 
                         levels = c("ISOTYPE", "DANA"), 
                         labels = c("Isotype", "Anti-LT\u03B13"))

data$sample <- factor(data$sample, 
                       levels = c("WT_LN_ISOTYPE", "WT_LN_DANA", "TNF_LN_ISOTYPE", "TNF_LN_DANA"), 
                       labels = c("WT Isotype", "WT Anti-LT\u03B13", "Tnf\u0394ARE Isotype", "Tnf\u0394ARE Anti-LT\u03B13"))

saveRDS(data, "mln_lt_blocking.rds")

cell_markers <- FindAllMarkers(data, only.pos = TRUE)
write.csv(cell_markers, "mln_cluster_markers.csv", row.names = FALSE)

