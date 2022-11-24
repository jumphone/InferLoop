#  /home/toolkit/tools/R4.2.0/bin/R


library(Signac)
library(Seurat)
library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v86) #hg38
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)


#############################################################i
library(monocle)
library(cicero)


pbmc=readRDS('scATAC_pbmc_final.rds')

library(data.table)

TAG='r0'
LOOP=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(1))$V1
SSN=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(2:(ncol(pbmc)+1)))
SSN=as.matrix(SSN)
rownames(SSN)=LOOP
colnames(SSN)=colnames(pbmc)


FEATURE="chr12-76254800-76255834.And.chr12-76709277-76710882"

pbmc[['isloop']]=CreateAssayObject(data=SSN)
DefaultAssay(pbmc)='isloop'

pbmc$type=readRDS('scATAC_type.rds')
Idents(pbmc)=pbmc$type

pdf('fxx_p03_loop_vlnPlot.pdf',width=2.5,height=2)
VlnPlot(pbmc, features=c(FEATURE),pt.size=0)+NoLegend()
dev.off()


#all.genes <- rownames(pbmc)
#pbmc <- ScaleData(pbmc, features = all.genes)
#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- RunUMAP(pbmc, dims = 1:10)
#DimPlot(pbmc, reduction = "umap")



#############################
SSN[which(SSN>4)]=4

pbmc[['isloop']]=CreateAssayObject(data=SSN)
DefaultAssay(pbmc)='isloop'

FEATURE="chr12-76254800-76255834.And.chr12-76709277-76710882"


pdf('fxx_p03_loop.pdf',width=5,height=5)
SpatialFeaturePlot(pbmc, features = c(FEATURE),pt.size.factor=3, alpha=0.9, crop=F)
SpatialFeaturePlot(pbmc, features = c(FEATURE),pt.size.factor=3, alpha=c(1,1), crop=F,stroke=0.25,image.alpha=1)
dev.off()



#ILS=pbmc[['isloop']]@data[which(rownames(pbmc[['isloop']]@data)==FEATURE),]













