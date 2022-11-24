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



conns=readRDS(file='/home/disk/database/data/scHIC/GSE130711_analysis/scATAC_conns.rds')
conns_pos=conns[which(conns[,3]>0),]



CUT=100000
LOOP=conns_pos[which(rank(-conns_pos[,3], ties.method='random')<= CUT*2 ),]
LOOP[,1]=stringr::str_replace_all(LOOP[,1],'_','-')
LOOP[,2]=stringr::str_replace_all(LOOP[,2],'_','-')
write.table(LOOP,  file='scATAC_loop_100k.txt', row.names=F,col.names=F,quote=F,sep='\t')


CUT=200000
LOOP=conns_pos[which(rank(-conns_pos[,3], ties.method='random')<= CUT*2 ),]
LOOP[,1]=stringr::str_replace_all(LOOP[,1],'_','-')
LOOP[,2]=stringr::str_replace_all(LOOP[,2],'_','-')
write.table(LOOP,  file='scATAC_loop_200k.txt', row.names=F,col.names=F,quote=F,sep='\t')


CUT=300000
LOOP=conns_pos[which(rank(-conns_pos[,3], ties.method='random')<= CUT*2 ),]
LOOP[,1]=stringr::str_replace_all(LOOP[,1],'_','-')
LOOP[,2]=stringr::str_replace_all(LOOP[,2],'_','-')
write.table(LOOP,  file='scATAC_loop_300k.txt', row.names=F,col.names=F,quote=F,sep='\t')


CUT=400000
LOOP=conns_pos[which(rank(-conns_pos[,3], ties.method='random')<= CUT*2 ),]
LOOP[,1]=stringr::str_replace_all(LOOP[,1],'_','-')
LOOP[,2]=stringr::str_replace_all(LOOP[,2],'_','-')
write.table(LOOP,  file='scATAC_loop_400k.txt', row.names=F,col.names=F,quote=F,sep='\t')




















