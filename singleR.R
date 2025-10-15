library(Seurat)
library(tidyverse)
library(SingleR)
library(celldex)
library(pheatmap)
library(ggplot2)
install.packages("viridis")
library(viridis)
library(tidyverse)
hd5_obj <- Read10X_h5(filename='20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5', use.names= TRUE, unique.features = TRUE)
pbmc_seurat <- CreateSeuratObject(counts = hd5_obj)
#QC and filtering
pbmc_seurat[['mt.percent']] <- PercentageFeatureSet(pbmc_seurat, pattern='^MT-')
View(pbmc_seurat@meta.data)
pbmc_seurat <- subset(pbmc_seurat, subset= nCount_RNA > 800 & nFeature_RNA > 500 & mt.percent < 10)
pbmc_seurat

#normalising, PCA, umap: not required for singleR, but for visualisation later
pbmc_seurat <- NormalizeData(pbmc_seurat)
pbmc_seurat <- FindVariableFeatures(pbmc_seurat)
pbmc_seurat <- ScaleData(pbmc_seurat)
pbmc_seurat <- RunPCA(pbmc_seurat)
pbmc_seurat <- FindNeighbors(pbmc_seurat, dims = 1:20)
pbmc_seurat <- FindClusters(pbmc_seurat)
pbmc_seurat <- RunUMAP(pbmc_seurat, dims = 1:20)
#get reference data from celldex
ref <- celldex::HumanPrimaryCellAtlasData()
#get the counts Data
pbmc_counts <- GetAssayData(pbmc_seurat, slot='counts')
#run singleR
pred <- SingleR(test = pbmc_counts,
        ref = ref,
        labels= ref$label.main)

#transfer the labels to the seurat object
pbmc_seurat$singleR.labels <- pred$labels[match(rownames(pbmc_seurat@meta.data), rownames(pred))]

#visualise the plot
DimPlot(pbmc_seurat, reduction='umap', group.by = 'singleR.labels')

#Annotation diagnostics
#based on scores- for a cell type label, the score must be clearly higher
plotScoreHeatmap(pred)
#based on delta values, lesser the delta value, less sure it is
plotDeltaDistribution(pred)
#based on comparison to unsupervised clustering
#create a table of each label, and how many cells are in each cluster- it does not know which label is right
tab <- table(Assigned=pred$labels, Clusters=pbmc_seurat@meta.data$seurat_clusters)
pheatmap(log10(tab+10), color=colorRampPalette(c('white','blue'))(10))
