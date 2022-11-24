# /home/toolkit/tools/R4.2.0/bin/R

setwd('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/tracks/')

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
#DimPlot(pbmc)

CLST=readRDS(file='/home/disk/database/data/scHIC/GSE130711_analysis/scATAC_clst.rds')

pbmc$clst=CLST

SSN=read.table('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/this.out.r0', sep='\t', row.names=1, header=TRUE,check.names=FALSE)



isLoop=matrix(0,nrow=nrow(SSN),ncol=ncol(pbmc))
rownames(isLoop)=rownames(SSN)
colnames(isLoop)=colnames(pbmc)

UC=unique(CLST)
i=1
while(i<=length(UC)){
    this_clst=UC[i]
    this_index1=which(CLST==this_clst)
    this_index2=which(colnames(SSN)==this_clst)
    isLoop[,this_index1]=SSN[,this_index2]
    print(i)
    i=i+1
    }

pbmc[['isloop']]=CreateAssayObject(data=isLoop)
DefaultAssay(pbmc)='isloop'

Idents(pbmc)=pbmc$type

DefaultAssay(pbmc)='isloop'
pbmc.markers <- FindAllMarkers(pbmc, slot='data',only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2, test.use = 'wilcox')

saveRDS(pbmc.markers, file='scATAC_isloop.markers.rds')


source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
REF=.generate_mean(isLoop, as.character(pbmc$type))


DREF=REF
DREF[which(DREF<0)]=0

i=1
while(i<=ncol(REF)){
    this_type=colnames(REF)[i]
    this_index=which( pbmc.markers$cluster == this_type)
    this_diff=pbmc.markers[this_index,]
    N=10000
    if(nrow(this_diff)<N){
        this_gene=this_diff$gene
        }else{
        this_gene=this_diff$gene[1:N]
        }
    #this_diff$gene[which(this_diff$p_val_adj<0.05)]
    DREF[which(!rownames(DREF) %in% this_gene),i]=0
    print(i)
    i=i+1
    }

length(which(DREF[,2]>0))


write.table(DREF,  file='scATAC_RefLoop.txt', row.names=T,col.names=T,quote=F,sep='\t')



































