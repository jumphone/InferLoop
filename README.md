<img src="https://fzhang.bioinfo-lab.com//img/inferloop_logo.jpg" width="250">

**InferLoop: leveraging single-cell chromatin accessibility for the signal of chromatin loop**

This tool is designed for inferring the loop signals of cell clusters (bins) in scATAC-seq data


# Demo data:

Key words: 10x genomics, scATAC-seq, PBMC

* [atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5)

* [atac_v1_pbmc_10k_singlecell.csv](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv)

* [atac_v1_pbmc_10k_fragments.tsv.gz](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz)

* [atac_v1_pbmc_10k_fragments.tsv.gz.tbi](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi)

* [pbmc_10k_v3.rds](https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds)

# Demo code:

[The official website of Signac](https://stuartlab.org/signac/)

## Section I, Using Signac to process the scATAC-seq data

    # /home/toolkit/tools/R4.2.0/bin/R
    setwd('/home/disk/database/data/ICS_scHIC/fzz_website_demo')

    library(Signac)
    library(Seurat)
    library(GenomeInfoDb)
    library(EnsDb.Hsapiens.v75)
    library(ggplot2)
    library(patchwork)
    set.seed(1234)
    
    counts <- Read10X_h5(filename = "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
    
    metadata <- read.csv(
        file = "atac_v1_pbmc_10k_singlecell.csv",
        header = TRUE,
        row.names = 1
        )
    
    chrom_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c(":", "-"),
        genome = 'hg19',
        fragments = 'atac_v1_pbmc_10k_fragments.tsv.gz',
        min.cells = 10,
        min.features = 200
        )

    pbmc <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "peaks",
        meta.data = metadata
        )
    
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
    annotations@seqnames@values=paste0('chr',annotations@seqnames@values)
    annotations@seqinfo@seqnames=paste0('chr',annotations@seqinfo@seqnames)
    annotations@seqinfo@genome=rep('hg19',length(annotations@seqinfo@genome))
    annotations@seqnames@values[which(annotations@seqnames@values=='chrMT')]='chrM'
    annotations@seqinfo@seqnames[which(annotations@seqinfo@seqnames=='chrMT')]='chrM'
    annotations@seqnames@values=factor(annotations@seqnames@values,levels=annotations@seqinfo@seqnames)
    
    Annotation(pbmc) <- annotations
    
    pbmc <- NucleosomeSignal(object = pbmc)
    pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
    pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
    pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
    
    pbmc <- subset(
        x = pbmc,
        subset = peak_region_fragments > 3000 &
        peak_region_fragments < 20000 &
        pct_reads_in_peaks > 15 &
        blacklist_ratio < 0.05 &
        nucleosome_signal < 4 &
        TSS.enrichment > 2
        )
    
    # Call peaks using MACS2
    
    peaks <- CallPeaks(pbmc, macs2.path = "/home/toolkit/local/bin/macs2")
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg19, invert = TRUE)
    
    saveRDS(peaks, 'peaks_macs2.rds')
    
    macs2_counts <- FeatureMatrix(
        fragments = Fragments(pbmc),
        features = peaks,
        cells = colnames(pbmc)
        )
    
    pbmc[["macs2"]] <- CreateChromatinAssay(
        counts = macs2_counts,
        fragments = 'atac_v1_pbmc_10k_fragments.tsv.gz',
        annotation = annotations
        )
    
    DefaultAssay(pbmc)='macs2'
    pbmc <- RunTFIDF(pbmc)
    pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
    pbmc <- RunSVD(pbmc)
    
    pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
    pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
    pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
    
    DimPlot(object = pbmc, label = TRUE) + NoLegend()

    gene.activities <- GeneActivity(pbmc)

    pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
    DefaultAssay(pbmc)='RNA'
    pbmc <- RunTFIDF(pbmc)
    
    pbmc_rna <- readRDS("pbmc_10k_v3.rds")

    transfer.anchors <- FindTransferAnchors(
        reference = pbmc_rna,
        query = pbmc,
        reduction = 'cca'
        )

     predicted.labels <- TransferData(
         anchorset = transfer.anchors,
         refdata = pbmc_rna$celltype,
         weight.reduction = pbmc[['lsi']],
         dims = 2:30
         )

    pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

    plot1 <- DimPlot(
        object = pbmc_rna,
        group.by = 'celltype',
        label = TRUE,
        repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

    plot2 <- DimPlot(
        object = pbmc,
        group.by = 'predicted.id',
        label = TRUE,
        repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

    plot1 + plot2
    
    saveRDS(pbmc, file='pbmc_signac.rds')
    
    

## Section II, Using Cicero to predict global loops

[The official website of Cicero](https://cole-trapnell-lab.github.io/cicero-release/docs/)  
    
    library(monocle)
    library(cicero)
    source('https://gitee.com/jumphone/public/raw/master/InferLoop.R')
    
    pbmc=readRDS(file='pbmc_signac.rds')
    DefaultAssay(pbmc)='macs2'
    
    indata=as.matrix(pbmc[['macs2']]@data)
    used_coords=pbmc@reductions$umap@cell.embeddings
    genome.df=inferloop.getGenomeDF.hg19()
    
    conns=inferloop.cicero(indata, used_coords, genome.df)
    
    saveRDS(conns, file='conns_cicero.rds')

## Section III, Using InferLoop.R to prepare input files
    
    source('https://gitee.com/jumphone/public/raw/master/InferLoop.R')
    
    pbmc=readRDS(file='pbmc_signac.rds')
    DefaultAssay(pbmc)='macs2'
    conns=readRDS(file='conns_cicero.rds')
    
    # The output loops of Cicero are not unique, and use top 400,000 loops will get top 200,000 unique loops
    inferloop.writeNet(conns, "net.txt", cut=400000) 
    
    indata=as.matrix(pbmc[['macs2']]@data)
    used_coords=pbmc@reductions$umap@cell.embeddings
    
    BIN=inferloop.generateBin(indata,used_coords, n=100)
    saveRDS(BIN, file='BIN.rds')
    
    CLST=BIN$clst
    write.table(BIN$mat,file='mat.txt', row.names=T,col.names=T,quote=F,sep='\t')
    
    
## Section IV, Using InferLoop to infer loop signals ( Python3 )
    
    mkdir output
    python3 inferloop/step0_uniqNet.py net.txt output/net_uniq.txt
    python3 inferloop/step1_buildIndex.py output/net_uniq.txt mat.txt output/mat.index
    python3 inferloop/step2_runInferLoop.py output/mat.index output/signal_mat.txt
    
    
## Section V, Generating inferred loop signals (ILS) of each cell type ( R )
    
    source('https://gitee.com/jumphone/public/raw/master/InferLoop.R')
    
    pbmc=readRDS(file='pbmc_signac.rds')
    DefaultAssay(pbmc)='macs2'
    
    MAT=inferloop.loadSignal('output/signal_mat.txt')
    
    CLST=readRDS('BIN.rds')$clst    
    MAT_CELL=inferloop.bin2cell(MAT, CLST)
    TYPE=pbmc$predicted.id
    MAT_TYPE=.generate_mean(MAT_CELL, TYPE)
    
    saveRDS(MAT_TYPE, 'signal_mat_type.rds')
    
    library(trackViewer)
    
    
    
    pbmc[['ILS']]=CreateAssayObject(data = MAT_CELL)
    
    
    
    

## Section VI, 
    
    

    
    
    
    
    
    
    
    
    
    
    
