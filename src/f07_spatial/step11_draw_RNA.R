


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
GLS=readRDS(paste0('./RNA/GLS_',TAG,'.rds'))




pbmc[['GLS']]=CreateAssayObject(data=GLS)
DefaultAssay(pbmc)='GLS'



pdf('fxx_p03_loop.pdf',width=5,height=5)
SpatialFeaturePlot(pbmc, features = c("Sptb"),pt.size.factor=3, alpha=0.8, crop=F)
dev.off()


Sptb_GLS=pbmc[['GLS']]@data[which(rownames(pbmc[['GLS']]@data)=='Sptb'),]
Sptb_GA=pbmc[['RNA']]@data[which(rownames(pbmc[['RNA']]@data)=='Sptb'),]

boxplot(Sptb_GLS~pbmc$type)
boxplot(Sptb_GA~pbmc$type)

t.test(Sptb_GLS[which(pbmc$type==5)],Sptb_GLS[which(pbmc$type!=5)])
t.test(Sptb_GA[which(pbmc$type==5)],Sptb_GA[which(pbmc$type!=5)])
























