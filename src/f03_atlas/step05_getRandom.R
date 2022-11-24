


#/home/toolkit/tools/R4.2.0/bin/R



library(Seurat)
library(stringr)
library(data.table)

LIST=read.table('LIST.txt',sep='\t')

pbmc=readRDS('./eQTL/scAtlas_pbmc.rds')
#ORGAN=t(matrix(unlist(str_split(colnames(pbmc),'_')),nrow=2))[,1]
#TYPE=colnames(pbmc)#t(matrix(unlist(str_split(colnames(pbmc),'_')),nrow=2))[,2]
#LOOP=read.table('./eQTL/SNP_TSS_500k.bed.loop.txt')

TAG='r0'
LOOP=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(1))$V1
SSN=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(2:(ncol(pbmc)+1)))
SSN=as.matrix(SSN)
rownames(SSN)=LOOP
colnames(SSN)=colnames(pbmc)

j=101
while(j<=500){
    set.seed(j)
    NROW=nrow(SSN)
    RSSN=SSN
    i=1
    while(i<=ncol(SSN)){
        RSSN[,i]=rnorm(n=NROW)
        i=i+1}
    rownames(RSSN)=rownames(SSN)
    write.table(RSSN, file=paste0('./ics_out/this.out.random.',j),sep='\t',quote=F,row.names=T,col.names=T)
    print(j)
    j=j+1
    }


