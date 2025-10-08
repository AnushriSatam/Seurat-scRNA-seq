library(Seurat)
library(SeuratDisk)
library(gridExtra)
library(tidyr)
library(ggplot2)
#get the location(subdirectory within the folder)
dirs <- list.dirs(path = "E:\\scRNA-seq pipeline\\integration", recursive = FALSE, full.names = FALSE)
#for loop to take out count matrix from each folder to make seurat objects
for(x in dirs)
{
  name <- gsub('_filtered_feature_bc_matrix', '',x)
  cts <-ReadMtx(mtx= paste0('E:\\scRNA-seq pipeline\\integration\\',x,'\\matrix.mtx.gz'),
        features=paste0('E:\\scRNA-seq pipeline\\integration\\',x,'\\features.tsv.gz'),
        cells=paste0('E:\\scRNA-seq pipeline\\integration\\',x,'\\barcodes.tsv.gz'))
  assign(name, CreateSeuratObject(counts=cts))
}
ls()
#combine the seurat objects
Merged_seurat  <- merge(HB17_background,y=c(HB17_PDX,HB17_tumor,HB30_PDX,HB30_tumor,HB53_background,HB53_tumor),
      add.cell.ids=ls()[3:9],
      project='HB')

 #QC and filtering
View(Merged_seurat@meta.data)


#Create sample column
Merged_seurat$sample <-rownames(Merged_seurat@meta.data)

#split the sample column
Merged_seurat@meta.data <- separate(Merged_seurat@meta.data,col=sample, into=c('Patient','Type','Barcode',sep="_"))

#calculate mitochondrial DNA
Merged_seurat[['percent.mt']] <- PercentageFeatureSet(Merged_seurat, pattern='^MT-')
FeatureScatter(nsclc_seurat_obj,feature1 = "nCount_RNA",feature2="nFeature_RNA")+
  geom_smooth(method="lm")
FeatureScatter(Merged_seurat,feature1 = "nCount_RNA",feature2="nFeature_RNA")+
  geom_smooth(method="lm")

#Filter out 
Merged_seurat_filtered <- subset(Merged_seurat,subset=nCount_RNA > 800 & nFeature_RNA > 500 & percent.mt < 10)

#when normalisation, variable features, PCA, uMAP was done, batch features due to patients was found so integration can be done:
#split the seurat object on basis of pateints
obj_list <- SplitObject(Merged_seurat_filtered, split.by='Patient')
for (i in 1:length(obj_list)){
  obj_list[[i]] <- NormalizeData(obj_list[[i]])
  obj_list[[i]] <- FindVariableFeatures(obj_list[[i]])
}
#select integration features
features <- SelectIntegrationFeatures(object.list = obj_list)

#find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, verbose = TRUE)

#integartion
seurat_integrated <- IntegrateData(anchorset=anchors)

#scaling
seurat_integrated <-ScaleData(seurat_integrated)

# PCA
seurat_integrated <- RunPCA(seurat_integrated)

#umap
seurat_integrated <- RunUMAP(seurat_integrated, dims= 1:15)

#plot
P3= DimPlot(seurat_integrated, reduction='umap', group.by='Patient')
P4= DimPlot(seurat_integrated, reduction='umap', group.by='Type')


grid.arrange(P3,P4,ncol=2)
