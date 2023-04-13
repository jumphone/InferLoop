

    # /home/toolkit/tools/R-4.2.0/bin/R
    
    setwd('/home/database/data/atacDevo/HiChIP_Cohesin_D0D7/ATAC/OUT')
    
    library(Signac)
    library(Seurat)
    library(GenomeInfoDb)
    library(EnsDb.Hsapiens.v75)
    library(ggplot2)
    library(patchwork)
    set.seed(1234)
    library(monocle)
    library(cicero)
    source('https://gitee.com/jumphone/public/raw/master/InferLoop.R')
    
    count_mat = inferloop.loadSignal('RMAT.txt')
    
    pbmc <- CreateSeuratObject(
        counts = count_mat,
        assay = "peaks"
        )
        
    pbmc <- RunTFIDF(pbmc)
    pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
    pbmc <- RunSVD(pbmc)
    pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
  
    DimPlot( object = pbmc, group.by='orig.ident')
    saveRDS(pbmc, file='pbmc_signac.rds')
    
    
    #####################################
    #####################################
    library(monocle)
    library(cicero)
    source('https://gitee.com/jumphone/public/raw/master/InferLoop.R')

    pbmc=readRDS(file='pbmc_signac.rds')
    DefaultAssay(pbmc)='peaks'

    indata=as.matrix(pbmc[['peaks']]@data)
    used_coords=pbmc@reductions$umap@cell.embeddings
    genome.df=inferloop.getGenomeDF.hg19()

    conns=inferloop.cicero(indata, used_coords, genome.df)

    saveRDS(conns, file='conns_cicero.rds')
    

    ######################################
    ######################################
    source('https://gitee.com/jumphone/public/raw/master/InferLoop.R')

    pbmc=readRDS(file='pbmc_signac.rds')
    DefaultAssay(pbmc)='peaks'
    conns=readRDS(file='conns_cicero.rds')

    # The output loops of Cicero are not unique, and use top 400,000 loops will get top 200,000 unique loops
    inferloop.writeNet(conns, "net.txt", cut=400000) 

    indata=as.matrix(pbmc[['peaks']]@data)
    used_coords=pbmc@reductions$umap@cell.embeddings

    BIN=inferloop.generateBin(indata,used_coords, n=100)
    saveRDS(BIN, file='BIN.rds')

    CLST=BIN$clst
    write.table(BIN$mat,file='mat.txt', row.names=T,col.names=T,quote=F,sep='\t')
   
    
    ######################################
    ######################################
    source('https://gitee.com/jumphone/public/raw/master/InferLoop.R')
    mat=inferloop.loadSignal('mat.txt')
    net=read.table('net.txt',sep='\t',header=F,row.names=NULL)
    net_uniq=inferloop.getUniqLoop(net)
    MAT=inferloop.inferLoopSignal(mat, net_uniq, r=0)
    write.table(MAT,file='signal_mat.txt', row.names=T,col.names=T,quote=F,sep='\t')
    
    
    ######################################
    ######################################
    source('https://gitee.com/jumphone/public/raw/master/InferLoop.R')
    library(Seurat)

    pbmc=readRDS(file='pbmc_signac.rds')
    DefaultAssay(pbmc)='peaks'

    MAT=inferloop.loadSignal('signal_mat.txt')

    CLST=readRDS('BIN.rds')$clst    
    MAT_CELL=inferloop.bin2cell(MAT, CLST)
    TYPE=pbmc$orig.ident
    MAT_TYPE=.generate_mean(MAT_CELL, TYPE)

    saveRDS(MAT_TYPE, 'signal_mat_type.rds')
    
    
   

