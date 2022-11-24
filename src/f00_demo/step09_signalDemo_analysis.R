# /home/toolkit/tools/R4.2.0/bin/R

library(Signac)
library(GenomicRanges)
library(Seurat)
#args <- commandArgs(TRUE)
#opath=args[1]
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')


pbmc=readRDS('/home/disk/database/data/scHIC/GSE130711_analysis/scATAC_pbmc.rds')
pbmc$type[which(pbmc$seurat_clusters %in% c(9))]='Unclear'

DISTAL=read.table('PDGFRA_TSS_hg19.1k.LOOP.distal.bed',sep='\t')
TSS=read.table('PDGFRA_TSS_hg19.1k.LOOP.tss.bed',sep='\t')

DefaultAssay(pbmc)='macs2'
pbmc <- RunTFIDF(pbmc)
DATA=pbmc[['macs2']]@data


UMAP=pbmc@reductions$umap@cell.embeddings
set.seed(123)
CLST=kmeans(UMAP,centers=200)$cluster
#CLST=readRDS(file='/home/disk/database/data/scHIC/GSE130711_analysis/scATAC_clst.rds')
TYPE=pbmc$type

mat=.generate_mean(DATA,CLST)
mat=mat[,order(as.numeric(colnames(mat)))]

tmp=table(TYPE, CLST)
type=rownames(tmp)[apply(-tmp,2,order)[1,]]



t_tag=paste0(TSS[,1],'-',TSS[,2],'-',TSS[,3])
d_tag=paste0(DISTAL[,1],'-',DISTAL[,2],'-',DISTAL[,3])
t_data=as.numeric(mat[which(rownames(mat) %in% t_tag),])
d_mat=mat[which(rownames(mat) %in% d_tag),]
d_mat_rsum=rowSums(d_mat)
#chr4-55306011-55307162
d_data=as.numeric(mat[which(rownames(mat) %in% names(which(d_mat_rsum==max(d_mat_rsum)))),])


pdf('f00p03_demo_pdgfra.pdf',width=3.5,height=4)

this_col=rep('grey60',length(t_data))
this_col[which(type=='OPC')]='red'
plot(t_data,d_data,col=this_col,pch=15,cex=2)
points(t_data[which(type=='OPC')],d_data[which(type=='OPC')],pch=15,col='red',cex=2)
dev.off()



a=t_data
b=d_data

out=rbind(a,b)
rownames(out)=c('a','b')
colnames(out)=c(1:ncol(out))
write.table(out,  file='signal_mat.txt', row.names=T,col.names=T,quote=F,sep='\t')

network=matrix(NA, nrow=1,ncol=2)
network[1,1]='a'
network[1,2]='b'
write.table(network,  file='signal_net.txt', row.names=F,col.names=F,quote=F,sep='\t')






saveRDS(type,file='signal_type.rds')
saveRDS(CLST,file='signal_CLST.rds')












