
# /home/toolkit/tools/R4.2.0/bin/R
library(Signac)
library(GenomicRanges)
library(Seurat)
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')



pbmc=readRDS('./scATAC_pbmc_final.rds')
#pbmc.isloop$type[which(pbmc.isloop$seurat_clusters %in% c(9))]='Unclear'
pbmc$type=readRDS('./scATAC_type.rds')

library(data.table)
    
    TAG='r0'
    LOOP=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(1))$V1
    SSN=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(2:(ncol(pbmc)+1)))
    SSN=as.matrix(SSN)
    rownames(SSN)=LOOP
    colnames(SSN)=colnames(pbmc)



pred.ref=SSN

LOOP2GENE=read.table("./RNA/LoopAndTSS.bed",header=F,sep='\t',row.names=NULL,check.names=F)
LOOP=paste0(LOOP2GENE[,1],'-',LOOP2GENE[,2],'-',LOOP2GENE[,3],'.And.',
            LOOP2GENE[,4],'-',LOOP2GENE[,5],'-',LOOP2GENE[,6])
GENE=paste0(LOOP2GENE[,11])

UGENE=unique(GENE)
OUT=matrix(0,ncol=ncol(pred.ref),nrow=length(UGENE))
rownames(OUT)=UGENE
colnames(OUT)=colnames(pred.ref)

i=1
while(i<=length(UGENE)){
    this_gene=UGENE[i]
    this_index=which(GENE==this_gene)
    this_loop=LOOP[this_index]
    this_loop_index=which(rownames(pred.ref) %in% this_loop)
    if(length(this_loop_index)>1){
        this_out=apply(pred.ref[this_loop_index,],2,sum)
        OUT[i,]=this_out
        }
    if(length(this_loop_index)==1){
        this_out=pred.ref[this_loop_index,]
        OUT[i,]=this_out
        }
    if(i %%500==1){print(i)}
    i=i+1}

saveRDS(OUT,paste0('./RNA/GLS_',TAG,'.rds'))





















