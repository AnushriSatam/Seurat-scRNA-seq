#load the library
library(Seurat)
library(tidyverse)
#load the dataset
nsclc_matrix <-Read10X_h5(filename='20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5')
str(nsclc_matrix) #it has 3 modalities
#choosing only the gene expression modality
cts <- nsclc_matrix$`Gene Expression`
#convert into seurat object
nsclc_seurat_obj <-CreateSeuratObject(counts=cts,project="NSCLC",min.cells=3,min.features = 200)
str(nsclc_seurat_obj)
view(nsclc_seurat_obj@meta.data) #nCount_RNA:total no of RNA present, nFeatureRNA: total no of unique genes present

#quality control to filter out low quality cells(those that have very low or very high (doublets) counts, high mitrochondrial gene content)
#calculate percentage of mt genes and add onto the meta data
nsclc_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(nsclc_seurat_obj,pattern="^MT-")
view(nsclc_seurat_obj@meta.data)

#voilin plot to visualise the total RNA, unique RNA, mt RNA in each cell, separately
VlnPlot(nsclc_seurat_obj,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

#look at it together (2 at a time)
FeatureScatter(nsclc_seurat_obj,feature1 = "nCount_RNA",feature2="nFeature_RNA")+
  geom_smooth(method="lm")

#filtering 
nsclc_seurat_obj <- subset(nsclc_seurat_obj, subset= nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalising
nsclc_seurat_obj <-NormalizeData(nsclc_seurat_obj)

#select genes that show high cell to cell variation
nsclc_seurat_obj <- FindVariableFeatures(nsclc_seurat_obj, selection.method = "vst",nfeatures=2000)

#top 10 highly variable genes
top_10 <- head(VariableFeatures(nsclc_seurat_obj),10)
top_10

#plot the variable genes and mark the top 10 ones
plot1 <- VariableFeaturePlot(nsclc_seurat_obj)
LabelPoints(plot=plot1,points=top_10, repel=TRUE)

#scaling 
all.genes <- rownames(nsclc_seurat_obj)
nsclc_seurat_obj <- ScaleData(nsclc_seurat_obj, features=all.genes)

#linear dimentionality reduction (PCA)
nsclc_seurat_obj <- RunPCA(nsclc_seurat_obj,features=VariableFeatures(object = nsclc_seurat_obj))

#determine dimentionality of the data(choose the PCAs that show high variance for clustering)
ElbowPlot(nsclc_seurat_obj)

#clustering

nsclc_seurat_obj <- FindNeighbors(nsclc_seurat_obj,dims= 1:15)
nsclc_seurat_obj <- FindClusters(nsclc_seurat_obj, resolution = c(0.1,0.3,0.5,0.7,1))
view(nsclc_seurat_obj@meta.data)

#view the clusters
DimPlot(nsclc_seurat_obj, group.by="RNA_snn_res.0.5",label=TRUE)
Idents(nsclc_seurat_obj) <- "RNA_snn_res.0.5"
Idents(nsclc_seurat_obj)

install.packages("umap")
#umap to visualise the clusters efficiently
nsclc_seurat_obj <- RunUMAP(nsclc_seurat_obj,dims= 1:15)
DimPlot(nsclc_seurat_obj, reduction="umap")
table(Idents(nsclc_seurat_obj))
