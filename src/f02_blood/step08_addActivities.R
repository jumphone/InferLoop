
library(Signac)
library(GenomicRanges)
library(Seurat)
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
#args <- commandArgs(TRUE)
#opath=args[1]

pbmc=readRDS('/home/disk/database/data/scHIC/Blood/analysis/scATAC_pbmc.rds')
#pbmc.isloop$type[which(pbmc.isloop$seurat_clusters %in% c(9))]='Unclear'
pbmc$type=readRDS('/home/disk/database/data/scHIC/Blood/analysis/scATAC_type.rds')
CLST=readRDS(file='/home/disk/database/data/scHIC/Blood/analysis/scATAC_clst.rds')


DefaultAssay(pbmc)='ATAC'

#gene.activities <- GeneActivity(pbmc)
#saveRDS(gene.activities,file='RNA/gene.activities.rds')

gene.activities=readRDS('RNA/gene.activities.rds')

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['ACT']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(pbmc)='ACT'
pbmc <- RunTFIDF(pbmc)
#pbmc <- NormalizeData(
#  object = pbmc,
#  assay = 'ACT',
#  normalization.method = 'LogNormalize',
#  scale.factor = median(pbmc$nCount_ACT)
#  )




is.loop=pbmc[['ACT']]@data
#############
    is.loop.ref=.generate_mean(is.loop, pbmc$type)
    pred.ref=as.matrix(is.loop.ref)
    #pred.ref=pred.ref[,which(!colnames(pred.ref) %in% c('Unclear'))]
    pred.ref=pred.ref[,which(colnames(pred.ref) %in% c('T','B','Mono'))]
    pred.ref=pred.ref[,order(colnames(pred.ref))]


tmp=readRDS('./RNA/COM_RNA.ref_r0.rds')
pred.ref=pred.ref[which(rownames(pred.ref) %in% rownames(tmp)),]


    colnames(pred.ref)=paste0(colnames(pred.ref),'_scATAC')
    saveRDS(pred.ref, './RNA/COM_ACT.rds' )














































