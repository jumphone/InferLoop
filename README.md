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

## Section 1, Using Signac to process the scATAC-seq data

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
    
    macs2_counts <- FeatureMatrix(
        fragments = Fragments(pbmc),
        features = peaks,
        cells = colnames(pbmc)
        )
    
    pbmc[["macs2"]] <- CreateChromatinAssay(
        counts = macs2_counts,
        fragments = 'atac_v1_pbmc_10k_fragments.tsv.gz',
        annotation = Annotation(pbmc)
        )
    
    DefaultAssay(pbmc)='macs2'
    pbmc <- RunTFIDF(pbmc)
    pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
    pbmc <- RunSVD(pbmc)
    
    pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
    pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
    pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

    pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
    DefaultAssay(pbmc)='RNA'
    pbmc <- RunTFIDF(pbmc)
    
    
    
    
    
    
    
    
