


#/home/toolkit/tools/R4.2.0/bin/R

library(Seurat)
library(Signac)
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')


LIST=read.table('LIST.txt',sep='\t')

MAT=c()

i=1
while(i<=nrow(LIST)){
    this_name=LIST[i,1]
    this_file=paste0('/home/disk/database/data/scATAC_ATLAS/raw/FRAG/',this_name,'.format/mat_snp_tss.type_agg.rds')
    this_mat=readRDS(this_file)
    if( this_name=='GSM5589390_ovary_SM-IOBHR_rep1_fragments.bed.gz'){
         MAT=cbind(MAT,this_mat) 
         MAT=MAT[,c(1:(ncol(MAT)-1))]
        }else{
         MAT=cbind(MAT,this_mat)
        }
    print(i)
    i=i+1
    }

saveRDS(MAT,'./eQTL/scAtlas_MAT.rds')


MAT=readRDS('./eQTL/scAtlas_MAT.rds')

DS_MAT=MAT

CSUM=colSums(MAT)
summary(CSUM)

#library('scuttle')
#PP=100000/CSUM
#PP[which(PP>1)]=1
#set.seed(123)
#DS_MAT=downsampleMatrix(MAT, PP)
#saveRDS(DS_MAT,'./disease/scAtlas_DS_MAT.rds')



pbmc <- CreateSeuratObject(
  counts = DS_MAT,
  assay = "ATAC",
  min.cells = 1
 )


DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)

OUT=pbmc[['ATAC']]@data
write.table(OUT,  file='./eQTL/scAtlas_normalizedData.txt', row.names=T,col.names=T,quote=F,sep='\t')

saveRDS(pbmc,file='./eQTL/scAtlas_pbmc.rds')


summary(colSums(OUT))







