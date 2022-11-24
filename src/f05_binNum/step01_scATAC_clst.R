 /home/toolkit/tools/R4.2.0/bin/R



library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)

library(Signac)
library(GenomicRanges)

pbmc=readRDS('/home/disk/database/data/scHIC/GSE130711_analysis/scATAC_pbmc.rds')
pbmc$type[which(pbmc$seurat_clusters %in% c(9))]='Unclear'
UMAP=pbmc@reductions$umap@cell.embeddings

DimPlot(pbmc,label=T,group.by='type')+NoLegend()


source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')

DATA= pbmc[['macs2']]@data


#############################################################
N=50
set.seed(123)
KM=kmeans(UMAP,centers=N)
CLST=KM$cluster
OUT=.generate_mean(DATA, CLST)
write.table(OUT,  file=paste0('./clst/scATAC_N',N,'_normalized_exp_clustered.txt'), row.names=T,col.names=T,quote=F,sep='\t')
saveRDS(CLST, file=paste0('./clst/scATAC_N',N,'_clst.rds'))


#############################################################
N=100
set.seed(123)
KM=kmeans(UMAP,centers=N)
CLST=KM$cluster
OUT=.generate_mean(DATA, CLST)
write.table(OUT,  file=paste0('./clst/scATAC_N',N,'_normalized_exp_clustered.txt'), row.names=T,col.names=T,quote=F,sep='\t')
saveRDS(CLST, file=paste0('./clst/scATAC_N',N,'_clst.rds'))


#############################################################
N=150
set.seed(123)
KM=kmeans(UMAP,centers=N)
CLST=KM$cluster
OUT=.generate_mean(DATA, CLST)
write.table(OUT,  file=paste0('./clst/scATAC_N',N,'_normalized_exp_clustered.txt'), row.names=T,col.names=T,quote=F,sep='\t')
saveRDS(CLST, file=paste0('./clst/scATAC_N',N,'_clst.rds'))

#############################################################
N=200
set.seed(123)
KM=kmeans(UMAP,centers=N)
CLST=KM$cluster
OUT=.generate_mean(DATA, CLST)
write.table(OUT,  file=paste0('./clst/scATAC_N',N,'_normalized_exp_clustered.txt'), row.names=T,col.names=T,quote=F,sep='\t')
saveRDS(CLST, file=paste0('./clst/scATAC_N',N,'_clst.rds'))















