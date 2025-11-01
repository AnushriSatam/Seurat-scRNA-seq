#calling the libraries
library(Seurat)
library(SeuratDisk)
library(harmony)
library(SeuratData)
library(tidyverse)
library(ggplot2)
library(umap)
#install dataset 
AvailableData()
InstallData("ifnb")
#load the data
ifnb <- LoadData("ifnb")
str(ifnb)
#QC
ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = '^MT-')
View(ifnb@meta.data)
FeatureScatter(ifnb, feature1 ="nFeature_RNA", feature2="nCount_RNA")+
  geom_smooth(method="lm")
#filtering
ifnb_filtered <- subset(ifnb,subset=nCount_RNA>800 & nFeature_RNA>200 & mito.percent<5)
#normalise
ifnb_filtered <- NormalizeData(ifnb_filtered)
ifnb_filtered <- FindVariableFeatures(ifnb_filtered)
ifnb_filtered <- ScaleData(ifnb_filtered)
#run PCA
ifnb_filtered <-RunPCA(ifnb_filtered)
ElbowPlot(ifnb_filtered)
#run umap
ifnb_filtered <- RunUMAP(ifnb_filtered, dims=1:15, reduction='pca')
DimPlot(ifnb_filtered, reduction = 'umap', group.by = 'stim')

#batch effect is seen, integration done by harmony
ifnb_harmony <- ifnb_filtered %>%
 RunHarmony(group.by.vars='stim',plot_convergence=FALSE)
ifnb_harmony@reductions
ifnb_harmony_embedding <- Embeddings (ifnb_harmony,"harmony")

#umap and clustering
ifnb_harmony <- ifnb_harmony %>%
    FindNeighbors(reduction='harmony',dims=1:15) %>%
    FindClusters(resolution=0.5) %>%
    RunUMAP(reduction='harmony', dims=1:15) %>%

#Plot
DimPlot(ifnb_harmony,reduction='umap', group.by='stim')
