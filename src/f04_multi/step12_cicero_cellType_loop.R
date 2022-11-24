# /home/toolkit/tools/R4.2.0/bin/R



library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
set.seed(1234)


library(Signac)
library(GenomicRanges)




pbmc=readRDS('./scATAC_pbmc.rds')
TYPE=pbmc$type


#############################################################i
getCicero<-function(pbmc,UMAP){
#################################################################
pbmc=pbmc
UMAP=UMAP
library(monocle)
library(cicero)
DefaultAssay(pbmc)='peaks'
pbmc <- RunTFIDF(pbmc)
indata=pbmc[['peaks']]@data
dim(indata)
rownames(indata)=stringr::str_replace_all(rownames(indata),'-','_')
cellinfo=pbmc@meta.data
peakinfo=as.data.frame(t(matrix(unlist(stringr::str_split(rownames(indata),'_')),nrow=3)))
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name=rownames(indata)
rownames(peakinfo) <- peakinfo$site_name
fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
input_cds <-  suppressWarnings(newCellDataSet(indata,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily=VGAM::binomialff(),
                            lowerDetectionLimit=0))
used_coords=UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = used_coords, size_factor_normalize = FALSE)

library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
genome <- genome[1:24]
genome.df <- data.frame("chr" = names(genome), "length" = genome)

source('/home/disk/database/data/ICS_scHIC/f04_multi/step02_scATAC_cicero_source.R')
conns <- .run_cicero_new(cicero_cds, genome.df)
#conns <- run_cicero(cicero_cds, genome.df)
head(conns)
return(conns)
##########################################################
}

DefaultAssay(pbmc)='peaks'

this_type='Neuron.Ex'
this_pbmc=subset(pbmc,cells=colnames(pbmc)[which(pbmc$type==this_type)])
this_umap=this_pbmc@reductions$umap@cell.embeddings
this_conns=getCicero(this_pbmc, this_umap)
saveRDS(this_conns, file=paste0('./cicero/',this_type,'_loop.rds'))


this_type='Neuron.In'
this_pbmc=subset(pbmc,cells=colnames(pbmc)[which(pbmc$type==this_type)])
this_umap=this_pbmc@reductions$umap@cell.embeddings
this_conns=getCicero(this_pbmc, this_umap)
saveRDS(this_conns, file=paste0('./cicero/',this_type,'_loop.rds'))


this_type='Astro'
this_pbmc=subset(pbmc,cells=colnames(pbmc)[which(pbmc$type==this_type)])
this_umap=this_pbmc@reductions$umap@cell.embeddings
this_conns=getCicero(this_pbmc, this_umap)
saveRDS(this_conns, file=paste0('./cicero/',this_type,'_loop.rds'))

this_type='ODC'
this_pbmc=subset(pbmc,cells=colnames(pbmc)[which(pbmc$type==this_type)])
this_umap=this_pbmc@reductions$umap@cell.embeddings
this_conns=getCicero(this_pbmc, this_umap)
saveRDS(this_conns, file=paste0('./cicero/',this_type,'_loop.rds'))

#############################



LOOP=read.table('scATAC_loop.txt.uniq.bed',header=FALSE,row.names=NULL,sep='\t')
#ALL_CONNS=readRDS('/home/disk/database/data/scHIC/Blood/analysis/scATAC_conns.rds')
#CUT=200000
#conns=ALL_CONNS
#conns_pos=conns[which(conns[,3]>0),]
#LOOP=conns_pos[which(rank(-conns_pos[,3], ties.method='random')<= CUT*2 ),]
LOOP=paste0(LOOP[,1],'_',LOOP[,2],'_',LOOP[,3],'_',LOOP[,4],'_',LOOP[,5],'_',LOOP[,6])


source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')

this_type='Neuron.Ex'
this_conns=readRDS(file=paste0('./cicero/',this_type,'_loop.rds'))
this_loop=paste0(this_conns[,1],'_',this_conns[,2])
this_index=which(this_loop %in% LOOP)
this_coa=this_conns[this_index,c(3)]
names(this_coa)=this_loop[this_index]
COM=cbind(this_coa,this_coa)

this_type='Neuron.In'
this_conns=readRDS(file=paste0('./cicero/',this_type,'_loop.rds'))
this_loop=paste0(this_conns[,1],'_',this_conns[,2])
this_index=which(this_loop %in% LOOP)
this_coa=this_conns[this_index,c(3)]
names(this_coa)=this_loop[this_index]
COM=.simple_combine(COM, cbind(this_coa,this_coa))$combine

this_type='Astro'
this_conns=readRDS(file=paste0('./cicero/',this_type,'_loop.rds'))
this_loop=paste0(this_conns[,1],'_',this_conns[,2])
this_index=which(this_loop %in% LOOP)
this_coa=this_conns[this_index,c(3)]
names(this_coa)=this_loop[this_index]
COM=.simple_combine(COM, cbind(this_coa,this_coa))$combine


this_type='ODC'
this_conns=readRDS(file=paste0('./cicero/',this_type,'_loop.rds'))
this_loop=paste0(this_conns[,1],'_',this_conns[,2])
this_index=which(this_loop %in% LOOP)
this_coa=this_conns[this_index,c(3)]
names(this_coa)=this_loop[this_index]
COM=.simple_combine(COM, cbind(this_coa,this_coa))$combine



COM=COM[,c(1,3,5,7)]
COM[is.na(COM)]=0
colnames(COM)=c('Neuron.Ex','Neuron.In','Astro','ODC')


RNAME=rownames(COM)
TMP=stringr::str_replace_all(RNAME,'_','-')
TMP=stringr::str_replace_all(TMP,'-chr','.And.chr')

rownames(COM)=TMP
COM[is.na(COM)]=0
COM=COM[,order(colnames(COM))]
colnames(COM)=paste0(colnames(COM),'_scATAC')


saveRDS(COM, file='./ics_out/COM_Cicero.rds')

SCOM=t(apply(COM,1,scale))
rownames(SCOM)=rownames(COM)
colnames(SCOM)=colnames(COM)


saveRDS(SCOM, file='./ics_out/SCOM_Cicero.rds')






