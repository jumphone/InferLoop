
#  /home/toolkit/tools/R4.2.0/bin/R


library(Signac)
library(Seurat)
library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v86) #hg38
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)


#############################################################i
library(monocle)
library(cicero)


pbmc=readRDS('scATAC_pbmc_final.rds')


UMAP=pbmc@reductions$umap@cell.embeddings


DefaultAssay(pbmc)='Spatial_ATAC'


indata=pbmc[['Spatial_ATAC']]@data
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

conns <- run_cicero(cicero_cds, genome.df)
head(conns)

saveRDS(conns, file='scATAC_conns.rds')

conns=readRDS(file='scATAC_conns.rds')


conns_pos=conns[which(conns[,3]>0),]


CUT=200000

LOOP=conns_pos[which(rank(-conns_pos[,3], ties.method='random')<= CUT*2 ),]
LOOP[,1]=stringr::str_replace_all(LOOP[,1],'_','-')
LOOP[,2]=stringr::str_replace_all(LOOP[,2],'_','-')

write.table(LOOP,  file='scATAC_loop.txt', row.names=F,col.names=F,quote=F,sep='\t')


DATA= pbmc[['Spatial_ATAC']]@data
write.table(DATA,  file='scATAC_normalized_exp_clustered.txt', row.names=T,col.names=T,quote=F,sep='\t')















