<img src="https://fzhang.bioinfo-lab.com//img/inferloop_logo.jpg" width="250">

**InferLoop: leveraging single-cell chromatin accessibility for the signal of chromatin loop**

This tool is designed for inferring the loop signals of cell clusters (bins) in scATAC-seq data


# Demo data:

Key words: 10x genomics, scATAC-seq, PBMC

* [atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5)

* [atac_v1_pbmc_10k_singlecell.csv](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv)

* [atac_v1_pbmc_10k_fragments.tsv.gz](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz)

* [atac_v1_pbmc_10k_fragments.tsv.gz.tbi](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi)

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
    
    
    
    
    
    
