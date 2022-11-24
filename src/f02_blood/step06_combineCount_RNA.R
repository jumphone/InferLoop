# /home/toolkit/tools/R4.2.0/bin/R

library(Signac)
library(GenomicRanges)
library(Seurat)
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
#args <- commandArgs(TRUE)
#opath=args[1]

pbmc.isloop=readRDS('/home/disk/database/data/scHIC/Blood/analysis/scATAC_pbmc.rds')
#pbmc.isloop$type[which(pbmc.isloop$seurat_clusters %in% c(9))]='Unclear'
pbmc.isloop$type=readRDS('/home/disk/database/data/scHIC/Blood/analysis/scATAC_type.rds')
CLST=readRDS(file='/home/disk/database/data/scHIC/Blood/analysis/scATAC_clst.rds')


pbmc.isloop[['EXP']]=pbmc.isloop[['RNA']]
DefaultAssay(pbmc.isloop)='EXP'
#pbmc.isloop <- RunTFIDF(pbmc.isloop)
pbmc.isloop <- NormalizeData(pbmc.isloop, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc.isloop)
pbmc.isloop <- ScaleData(pbmc.isloop, features = all.genes)

is.loop=pbmc.isloop[['EXP']]@scale.data
#############
    is.loop.ref=.generate_mean(is.loop, pbmc.isloop$type)
    pred.ref=as.matrix(is.loop.ref)
    #pred.ref=pred.ref[,which(!colnames(pred.ref) %in% c('Unclear'))]
    pred.ref=pred.ref[,which(colnames(pred.ref) %in% c('T','B','Mono'))]
    pred.ref=pred.ref[,order(colnames(pred.ref))]
    colnames(pred.ref)=paste0(colnames(pred.ref),'_scRNA')
    saveRDS(pred.ref, './RNA/COM_RNA.rds' )


#plot(pred.ref[,1],pred.ref[,2])
#cor(pred.ref)





