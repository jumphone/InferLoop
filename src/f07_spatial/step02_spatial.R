
#  /home/toolkit/tools/R4.2.0/bin/R


library(Signac)
library(Seurat)
library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v86) #hg38
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

pbmc=readRDS('scATAC_pbmc.rds')

image=Read10X_Image('/home/disk/database/data/ICS_scHIC/fxx_spatial/download/GSM5238386_ME13_50um_spatial')
image <- image[Cells(x = pbmc)]
assay = "Spatial_RNA"
slice = "slice1"
data=pbmc[['RNA']]@counts
#object <- CreateSeuratObject(counts = data, assay = assay)
#DefaultAssay(object = image) <- assay
#object[[slice]] <- image
pbmc[[assay]]=CreateAssayObject(counts=data)
DefaultAssay(pbmc = image) <- assay
pbmc[[slice]] <- image


VlnPlot(
  object = pbmc,
  features = c('nCount_peaks','nFeature_peaks'),
  pt.size = 0.1,
  ncol = 2
)

pbmc=subset(pbmc, subset =nCount_peaks<200000 & nCount_peaks > 20000)
pbmc=subset(pbmc, cells=rownames(pbmc@images$slice1@coordinates))


#pbmc <- SCTransform(pbmc, assay = "Spatial_RNA", verbose = TRUE)
DefaultAssay(pbmc)='Spatial_RNA'
pbmc <- RunTFIDF(pbmc)
#SpatialFeaturePlot(pbmc, features = c("Sptb"),pt.size.factor=4, alpha=0.7)


image=Read10X_Image('/home/disk/database/data/ICS_scHIC/fxx_spatial/download/GSM5238386_ME13_50um_spatial')
image <- image[Cells(x = pbmc)]
assay = "Spatial_ATAC"
slice = "slice1"
data=pbmc[['macs2']]@counts
#object <- CreateSeuratObject(counts = data, assay = assay)
#DefaultAssay(object = image) <- assay
#object[[slice]] <- image
pbmc[[assay]]=CreateAssayObject(counts=data)
DefaultAssay(pbmc = image) <- assay
pbmc[[slice]] <- image

DefaultAssay(pbmc)='Spatial_ATAC'
pbmc <- RunTFIDF(pbmc)

saveRDS(pbmc,file='scATAC_pbmc_final.rds')













