

#/home/toolkit/tools/R4.2.0/bin/R

library(Seurat)
library(Signac)
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
library(data.table)
library(stringr)



pbmc=readRDS('./eQTL/scAtlas_pbmc.rds')
SNP_TSS=read.table('./eQTL/SNP_TSS_500k.bed.loop.txt',sep='\t',row.names=NULL,header=FALSE)
SNP_TSS_TAG=SNP_TSS[,3]
SNP_TSS_GENE=SNP_TSS[,4]
INFO=cbind(paste0(SNP_TSS[,1],'.And.',SNP_TSS[,2]),SNP_TSS_TAG,SNP_TSS_GENE)
rownames(INFO)=INFO[,1]

PAIR=paste0(INFO[,2],'_',INFO[,3])

MATCH=read.table('MATCH.txt',sep='\t')


####################
source('step06_analysis_source.R')
##########################
OUT=c()
##########################
TAG='mean'
this_out=.getSub(TAG)
OUT=cbind(OUT,this_out)
##########################
TAG='hmean'
this_out=.getSub(TAG)
OUT=cbind(OUT,this_out)
##########################
TAG='ori'
this_out=.getSub(TAG)
OUT=cbind(OUT,this_out)
##########################
TAG='cicero_scale'
this_out=.getSub(TAG)
OUT=cbind(OUT,this_out)
##########################
TAG='r50'
this_out=.getSub(TAG)
OUT=cbind(OUT,this_out)
##########################
TAG='r0'
this_out=.getSub(TAG)
OUT=cbind(OUT,this_out)



saveRDS(OUT,'compare_real_scaleCicero.rds')
#saveRDS(OUT,'compare_real.rds')












