
# /home/toolkit/tools/R4.2.0/bin/R
p
p
p
p
p
library(Signac)
library(GenomicRanges)
library(Seurat)
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')




.getRNA<-function(pred.ref){

#pred.ref=readRDS('./ics_out/COM_pred.ref_r0.rds')
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

return(OUT)
}

####################
pred.ref=readRDS('./ics_out/COM_pred.ref_mean.rds')
rna.ref=.getRNA(pred.ref)
saveRDS(rna.ref,'./RNA/COM_RNA.ref_mean.rds')
###############
pred.ref=readRDS('./ics_out/COM_pred.ref_hmean.rds')
rna.ref=.getRNA(pred.ref)
saveRDS(rna.ref,'./RNA/COM_RNA.ref_hmean.rds')
#################
pred.ref=readRDS('./ics_out/COM_pred.ref_ori.rds')
rna.ref=.getRNA(pred.ref)
saveRDS(rna.ref,'./RNA/COM_RNA.ref_ori.rds')
##############
pred.ref=readRDS('./ics_out/COM_pred.ref_r50.rds')
rna.ref=.getRNA(pred.ref)
saveRDS(rna.ref,'./RNA/COM_RNA.ref_r50.rds')
#############
pred.ref=readRDS('./ics_out/COM_pred.ref_r0.rds')
rna.ref=.getRNA(pred.ref)
saveRDS(rna.ref,'./RNA/COM_RNA.ref_r0.rds')

##############





















