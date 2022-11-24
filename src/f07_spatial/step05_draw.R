


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

library(data.table)

TAG='r0'
LOOP=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(1))$V1
SSN=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(2:(ncol(pbmc)+1)))
SSN=as.matrix(SSN)
rownames(SSN)=LOOP
colnames(SSN)=colnames(pbmc)

DefaultAssay(pbmc)='Spatial_RNA'

#SpatialFeaturePlot(pbmc, features = c("Sptb"),pt.size.factor=3, alpha=c(0,1), crop=F,stroke=0)

SpatialFeaturePlot(pbmc, features = c("Sptb"),pt.size.factor=3, alpha=c(0,1), crop=F,stroke=1,image.alpha=0.5)

pdf('fxx_p01_sptb.pdf',width=5,height=5)
SpatialFeaturePlot(pbmc, features = c("Sptb"),pt.size.factor=3, alpha=0.8, crop=F)
SpatialFeaturePlot(pbmc, features = c("Sptb"),pt.size.factor=3, alpha=c(1,1), crop=F,stroke=0.25,image.alpha=1)
dev.off()


UMAP=pbmc@reductions$umap@cell.embeddings
DefaultAssay(pbmc)='Spatial_ATAC'
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

pdf('fxx_p02_cluster.pdf',width=5,height=5)
SpatialDimPlot(pbmc, label = TRUE, label.size = 10, crop=F,pt.size.factor=3)+NoLegend()
dev.off()

saveRDS(pbmc$Spatial_ATAC_snn_res.0.8,file='scATAC_type.rds')



pbmc$type=readRDS('scATAC_type.rds')

pbmc[['isloop']]=CreateAssayObject(data=SSN)
DefaultAssay(pbmc)='isloop'
Idents(pbmc)=pbmc$type
pbmc.markers <- FindAllMarkers(pbmc, slot='data',only.pos = TRUE, min.pct = 0, logfc.threshold = 0, test.use = 't')


saveRDS(pbmc.markers, file='scATAC_isloop.markers.rds')

#########################################

pbmc.markers=readRDS( file='scATAC_isloop.markers.rds')


source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
REF=.generate_mean(SSN, as.character(pbmc$type))

DREF=REF
i=1
while(i<=ncol(REF)){
    this_type=colnames(REF)[i]
    this_index=which( pbmc.markers$cluster == this_type)
    this_diff=pbmc.markers[this_index,]
    N=10000
    if(nrow(this_diff)<N){
        this_gene=this_diff$gene
        }else{
        this_gene=this_diff$gene[1:N]
        }
    #this_diff$gene[which(this_diff$p_val_adj<0.05)]
    DREF[which(!rownames(DREF) %in% this_gene),i]=0
    print(i)
    i=i+1
    }

length(which(DREF[,2]>0))

write.table(DREF,  file='scATAC_RefLoop.txt', row.names=T,col.names=T,quote=F,sep='\t')




















