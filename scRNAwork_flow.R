rm(list = ls()) 
library(Seurat)
library(SingleR)
library(celldex)
library(harmony)
library(dplyr)
library(tidyverse)
library(patchwork)
load("seurat.Rdata")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
scobj <- subset(scobj, subset = nFeature_RNA > 500 & percent.rb < 50 & percent.mt < 20)
scobj <- subset(scobj, subset = nFeature_RNA < 5000 & nCount_RNA < 40000)

scobj <- SCTransform(scobj, vars.to.regress = c("percent.mt"), verbose = FALSE)
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj), reduction.name = "pca")
ElbowPlot(scobj, reduction = "pca", ndims = 50)
xx <- cumsum(scobj[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.9) 
ndim <- 34 

scobj <- RunHarmony(scobj, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
scobj <- FindNeighbors(scobj,reduction = 'harmony',dims = 1:34)
scobj <- RunUMAP(scobj,reduction = "harmony",dims = 1:34)
resolutions <- c(0.6,0.7, 0.8, 0.9,1)
for (res in resolutions) {
  scobj <- FindClusters(scobj, resolution = res, verbose = FALSE)
  cluster_col_name <- paste0("seurat_clusters_res_", res)
  scobj[[cluster_col_name]] <- Idents(scobj)
}

for (res in resolutions) {
  cluster_col_name <- paste0("seurat_clusters_res_", res)
  p <- DimPlot(scobj, reduction = "umap", group.by = cluster_col_name, label = TRUE) + 
    ggtitle(paste("Resolution:", res))
  print(p)
}
Idents(scobj)<-"seurat_clusters_res_0.6"
DimPlot(scobj,reduction = "umap",label = T,group.by = "seurat_clusters_res_0.6")
markers <- FindAllMarkers(object = scobj, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,file = "markergene.csv")

genes <- c("CD19","CD79A","CD79B","MS4A1",'TPSAB1','TPSB2','CD68','CD14',
           'CDH5','VWF','PDGFRB','ACTA2','MOG','MBP','CD3D','CD3E',"FGFBP2","NKG7",'PDGFRA','SOX2',"OLIG1","GFAP")

DotPlot(scobj, features = genes) + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

FeaturePlot(scobj,features = c("CD19","CD79A","CD79B","MS4A1"),reduction = "umap")
FeaturePlot(scobj,features = c('TPSAB1','TPSB2'),reduction = "umap")
FeaturePlot(scobj,features = c('CD68','CD14'),reduction = "umap")
FeaturePlot(scobj,features = c('CDH5','VWF'),reduction = "umap")
FeaturePlot(scobj,features = c('PDGFRB','ACTA2'),reduction = "umap")
FeaturePlot(scobj,features = c('MOG','MBP'),reduction = "umap")
FeaturePlot(scobj,features = c('CD3D','CD3E'),reduction = "umap")
FeaturePlot(scobj,features = c("FGFBP2","NKG7"),reduction = "umap")
FeaturePlot(scobj,features = c('PDGFRA','SOX2',"OLIG1","GFAP"),reduction = "umap")

cluster_renaming <- c(
  "24" = "B cell", 
  "26" = "Mast cell", 
  "16" = "Oligo", 
  "6" = "T cell", 
  "20" = "NK cell", 
  "15" = "Pericyte", 
  "18" = "Endo", 
  "9" = "Glioma", "23" = "Glioma", "8" = "Glioma", "11" = "Glioma", 
  "14" = "Glioma", "13" = "Glioma", "10" = "Glioma", "7" = "Glioma", 
  "2" = "Glioma","3" = "Glioma","21" = "Glioma","22" = "Glioma",
  "19" = "Myeloid", "12" = "Myeloid", "25" = "Myeloid", "5" = "Myeloid", 
  "17" = "Myeloid", "0" = "Myeloid", "1" = "Myeloid", "4" = "Myeloid")

scobj <- RenameIdents(scobj, cluster_renaming)

scobj$celltype <- Idents(scobj)

signature<-c("COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL6A2",
             "COL6A3","COL18A1","ACTA2","AEBP1","BGN","CD248","ITGA5","LAMB1","LOXL1",
             "LUM","MMP9","NID2","PCOLCE","SERPINH1","TAGLN","VIM","TGFBI")
signature <- intersect(signature, rownames(scobj))

scobj <- AddModuleScore(
  object = scobj,
  features = list(signature),
  name = "signature"
)