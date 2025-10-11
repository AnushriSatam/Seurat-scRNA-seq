library(Seurat)
library(tidyverse)
library(metap)
library(multtest)
#load the data-already integrated
ifnb_harmony <- readRDS('ifnb_harmony.rds') 
str(ifnb_harmony)
View(ifnb_harmony@meta.data)
#visualise the clusters
DimPlot(ifnb_harmony,reduction='umap',group.by = 'seurat_clusters',label = TRUE)
#check the default assay. in case of harmony, no serparate assay is made (no corrected expression matrix given) as in the case of seurat. so the default assay will be 'RNA', but in case of seurat it will be 'integrated'
#check the default assay
DefaultAssay(ifnb_harmony)
#convert the default assay to RNA if not already
#DefaultAssay(ifnb_harmony) <- 'RNA'
#identify the cell types within the clusters
#findallmarkers would be more suitable if there was only one condition. here we have 2 conditions under stim
FindAllMarkers(ifnb_harmony,
               logfc.threshold = 0.25,
               min.pct= 0.1,
               only.pos= TRUE,               #want only upregualted one
               test.use = 'DeSeq2',         #Deseq2 needs counts
               slot = 'counts')

#use findconservedmarkers- for a cluster (ident.1) against others or a specific one (if specific: ident.2). it must be used in case of multiple conditions, if not specific 
markers_3 <- FindConservedMarkers(ifnb_harmony,
                     ident.1 = 3,
                     grouping.var = 'stim')
#visualise top features
FeaturePlot(ifnb_harmony,features=c('FCGR3A'), min.cutoff='q10')

#rename the idents (as in the name of the cells, which is based on the cluster number to a cell name based on the  marker)
Idents(ifnb_harmony)
ifnb_harmony <- RenameIdents(ifnb_harmony, '3'='CD16 mono')

#Visualise
DimPlot(ifnb_harmony,reduction='umap',label=T)

#in this dataset, the annotations have been provided, which is not an ideal scenario.
#so we can easily set the ident names as the seurat_annotation
Idents(ifnb_harmony) <-ifnb_harmony@meta.data$seurat_annotations
DimPlot(ifnb_harmony, reduction='umap', label=T)

#you want the cell type with the conditon. concatenation:
ifnb_harmony$celltype.cnd <- paste0(ifnb_harmony$seurat_annotations,'_',ifnb_harmony$stim)
View(ifnb_harmony@meta.data)

#set the idents value to this new column
Idents(ifnb_harmony) <- ifnb_harmony$celltype.cnd
DimPlot(ifnb_harmony,reduction = 'umap',label=T)

#comapre the control and stimulated cells within the same cell type
b.interferon.respnse <- FindMarkers(ifnb_harmony, ident.1 = 'CD16 Mono_CTRL', ident.2 = 'CD16 Mono_STIM')
View(b.interferon.respnse)

#plot conserved features within cluster, and the differentially expressed ones withn the control and stimulated of the cluster

FeaturePlot(ifnb_harmony, features=c('FCGR3A','IFIT1'), split.by = 'stim', min.cutoff = 'q10')

#findconservedmarker tells what makes CD16 cells CD16

