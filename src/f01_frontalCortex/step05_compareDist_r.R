
# /home/toolkit/tools/R4.2.0/bin/R

library(Signac)
library(GenomicRanges)
library(Seurat)
#args <- commandArgs(TRUE)
#opath=args[1]
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')


.getYesNo<-function(DType,CType){
    DType=DType
    CType=CType
    i.yes=c()
    j.yes=c()
    i.no=c()
    j.no=c()
    i=1
    while(i<length(DType)){
        j=i+1
        while(j<=length(DType)){
            this_i_dtype=DType[i]
            this_j_dtype=DType[j]
            this_i_ctype=CType[i]
            this_j_ctype=CType[j]
            if(this_i_dtype != this_j_dtype & this_i_ctype == this_j_ctype){i.yes=c(i.yes, i);j.yes=c(j.yes, j)}
            if(this_i_ctype != this_j_ctype){i.no=c(i.no, i);j.no=c(j.no, j)}
            j=j+1}
        i=i+1}
    
    this_out=list()
    this_out$i.yes=i.yes
    this_out$j.yes=j.yes
    this_out$i.no=i.no
    this_out$j.no=j.no
    return(this_out)
    }


.getSub <-function(real.ref,pred.ref){
    COM=.simple_combine(real.ref,pred.ref)
    D1=COM$exp_sc_mat1
    D2=COM$exp_sc_mat2
    ALL=cbind(D1,D2)
    COR=cor(ALL,method='pearson')
    mat=COR
    DType=t(matrix(unlist(stringr::str_split(colnames(mat),'_')),nrow=2))[,2]
    CType=t(matrix(unlist(stringr::str_split(colnames(mat),'_')),nrow=2))[,1]
    this_out=.getYesNo(DType,CType)
    this_yes=diag(mat[this_out$i.yes, this_out$j.yes])
    this_no=diag(mat[this_out$i.no, this_out$j.no])
    #this_no_head=this_no[order(-this_no)][1:length(this_yes)]
    this_no_head=this_no[order(-this_no)][1]
    this_out=list()
    this_out$yes=this_yes
    this_out$no=this_no
    this_out$no_head=this_no_head
    #this_sub=mean(this_yes)-mean(this_no)
    return(this_out)
    }

#real.ref=readRDS('/home/disk/database/data/scHIC/GSE130711_analysis/NEW/COM_real.ref.rds')
real.ref=readRDS('./ics_out/COM_real.ref.rds')
real.ref=real.ref[which(rowSums(real.ref)>0),]

OUT=list()
SUB=c()
################################################################
pred.ref=readRDS('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/COM_pred.ref_nr200.rds')
this_out=.getSub(real.ref,pred.ref)
OUT[[1]]=this_out$yes
OUT[[2]]=this_out$no_head
################################################################
pred.ref=readRDS('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/COM_pred.ref_nr150.rds')
this_out=.getSub(real.ref,pred.ref)
OUT[[3]]=this_out$yes
OUT[[4]]=this_out$no_head
################################################################
pred.ref=readRDS('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/COM_pred.ref_nr100.rds')
this_out=.getSub(real.ref,pred.ref)
OUT[[5]]=this_out$yes
OUT[[6]]=this_out$no_head
################################################################
pred.ref=readRDS('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/COM_pred.ref_nr50.rds')
this_out=.getSub(real.ref,pred.ref)
OUT[[7]]=this_out$yes
OUT[[8]]=this_out$no_head
################################################################
pred.ref=readRDS('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/COM_pred.ref_r0.rds')
this_out=.getSub(real.ref,pred.ref)
OUT[[9]]=this_out$yes
OUT[[10]]=this_out$no_head
################################################################
pred.ref=readRDS('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/COM_pred.ref_r50.rds')
this_out=.getSub(real.ref,pred.ref)
OUT[[11]]=this_out$yes
OUT[[12]]=this_out$no_head
################################################################
pred.ref=readRDS('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/COM_pred.ref_r100.rds')
this_out=.getSub(real.ref,pred.ref)
OUT[[13]]=this_out$yes
OUT[[14]]=this_out$no_head
################################################################
pred.ref=readRDS('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/COM_pred.ref_r150.rds')
this_out=.getSub(real.ref,pred.ref)
OUT[[15]]=this_out$yes
OUT[[16]]=this_out$no_head
################################################################
pred.ref=readRDS('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/COM_pred.ref_r200.rds')
this_out=.getSub(real.ref,pred.ref)
OUT[[17]]=this_out$yes
OUT[[18]]=this_out$no_head

################################################################

#boxplot(OUT,col=c('indianred1','royalblue1'))

SUB=c()
SUB=c(SUB,mean(OUT[[1]])-mean(OUT[[2]]))
SUB=c(SUB,mean(OUT[[3]])-mean(OUT[[4]]))
SUB=c(SUB,mean(OUT[[5]])-mean(OUT[[6]]))
SUB=c(SUB,mean(OUT[[7]])-mean(OUT[[8]]))
SUB=c(SUB,mean(OUT[[9]])-mean(OUT[[10]]))
SUB=c(SUB,mean(OUT[[11]])-mean(OUT[[12]]))
SUB=c(SUB,mean(OUT[[13]])-mean(OUT[[14]]))
SUB=c(SUB,mean(OUT[[15]])-mean(OUT[[16]]))
SUB=c(SUB,mean(OUT[[17]])-mean(OUT[[18]]))


print(SUB)
setwd('/home/disk/database/data/ICS_scHIC/f01_frontalCortex')

pdf('p03_sub_barplot_cortex_r.pdf',width=10,height=5)

barplot(SUB,col=c(rep('grey60',4),'indianred1',rep('grey60',4)),
        lwd=1,border = NA, space = 2)
abline(h=0,lwd=1.5)
dev.off()



















