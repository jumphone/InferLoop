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


############################
get_pred <- function(TAG){    
    SSN=read.table(paste0("./ics_out/this.out.",TAG),header=T,sep='\t',row.names=1,check.names=F)
    isLoop=matrix(0,nrow=nrow(SSN),ncol=ncol(pbmc.isloop))
    rownames(isLoop)=rownames(SSN)
    colnames(isLoop)=colnames(pbmc.isloop)
    ##########
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
    is.loop=isLoop
    #############
    is.loop.ref=.generate_mean(is.loop, pbmc.isloop$type)
    pred.ref=as.matrix(is.loop.ref)
    #pred.ref=pred.ref[,which(!colnames(pred.ref) %in% c('Unclear'))]
    pred.ref=pred.ref[,which(colnames(pred.ref) %in% c('T','B','Mono'))]
    pred.ref=pred.ref[,order(colnames(pred.ref))]
    colnames(pred.ref)=paste0(colnames(pred.ref),'_scATAC')
    saveRDS(pred.ref, paste0('./ics_out/COM_pred.ref_',TAG,'.rds') )
    }
#############################


get_pred('mean')
get_pred('hmean')
get_pred('ori')
get_pred('r0')
get_pred('r50')









