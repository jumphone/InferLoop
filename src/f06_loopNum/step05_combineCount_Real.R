
# /home/toolkit/tools/R4.2.0/bin/R

library(Signac)
library(GenomicRanges)
library(Seurat)
#args <- commandArgs(TRUE)
#opath=args[1]

.single=function(TAG){

d1=read.table(paste0('./coolCount_',TAG,'/Astro.coolCount.txt'),sep='\t',header=FALSE,row.names=NULL)
d2=read.table(paste0('./coolCount_',TAG,'/Neuron.Ex.coolCount.txt'),sep='\t',header=FALSE,row.names=NULL)
d3=read.table(paste0('./coolCount_',TAG,'/Neuron.In.coolCount.txt'),sep='\t',header=FALSE,row.names=NULL)
d4=read.table(paste0('./coolCount_',TAG,'/MG.coolCount.txt'),sep='\t',header=FALSE,row.names=NULL)
d5=read.table(paste0('./coolCount_',TAG,'/ODC.coolCount.txt'),sep='\t',header=FALSE,row.names=NULL)
d6=read.table(paste0('./coolCount_',TAG,'/OPC.coolCount.txt'),sep='\t',header=FALSE,row.names=NULL)

DATA=cbind(d1[,7],d2[,7],d3[,7],d4[,7],d5[,7],d6[,7])
colnames(DATA)=c('Astro','Neuron.Ex','Neuron.In','MG','ODC','OPC')
RNAME=paste0(d1[,1],'-',d1[,2],'-',d1[,3],'.And.',d1[,4],'-',d1[,5],'-',d1[,6])
rownames(DATA)=RNAME

pbmc.ref <- CreateSeuratObject(
  counts = DATA,
  assay = "peaks"
  )
pbmc.ref <-  RunTFIDF(pbmc.ref)

NDATA=pbmc.ref[['peaks']]@data
real.ref=as.matrix(NDATA)
real.ref=real.ref[,which(!colnames(real.ref) %in% c('Unclear'))]
real.ref=real.ref[,order(colnames(real.ref))]
colnames(real.ref)=paste0(colnames(real.ref),'_scHIC')

saveRDS(real.ref,paste0('./coolCount_',TAG,'/COM_real.ref.rds'))
}

.single('100k')
.single('200k')
.single('300k')
.single('400k')






