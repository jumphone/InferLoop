
 /home/toolkit/tools/R4.2.0/bin/R


setwd('/home/disk/database/data/ICS_scHIC/f04_multi')

library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
#diyn.load('/usr/local/hdf5-1.8.17/lib/libhdf5_hl.so.10')

set.seed(1234)


rna <- Read10X("./data/rna/", gene.column = 1)
atac <- Read10X("./data/atac/", gene.column = 1)
fragments <- "/home/disk/database/data/ICS_scHIC/f04_multi/data/fragments.sort.bed.gz"

# create a Seurat object and add the assays
pbmc <- CreateSeuratObject(counts = rna)
pbmc[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = fragments
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)


annotations@seqnames@values=paste0('chr',annotations@seqnames@values)
annotations@seqinfo@seqnames=paste0('chr',annotations@seqinfo@seqnames)
annotations@seqinfo@genome=rep('mm10',length(annotations@seqinfo@genome))
annotations@seqnames@values[which(annotations@seqnames@values=='chrMT')]='chrM'
annotations@seqinfo@seqnames[which(annotations@seqinfo@seqnames=='chrMT')]='chrM'
annotations@seqnames@values=factor(annotations@seqnames@values,levels=annotations@seqinfo@seqnames)
#saveRDS(annotations,'/home/database/annotation/mm10/mm10_signac_ucsc_annotations.rds')


annotations=readRDS('/home/database/annotation/mm10/mm10_signac_ucsc_annotations.rds')

# add the gene information to the object
Annotation(pbmc[["ATAC"]]) <- annotations



DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)
pbmc$blacklist_fraction <- FractionCountsInRegion(
  object = pbmc,
  assay = 'ATAC',
  regions = blacklist_mm10
)

pbmc <- subset(
  x = pbmc,
  subset = blacklist_fraction < 0.03 &
    TSS.enrichment < 20 &
    nCount_RNA > 800 &
    nCount_ATAC > 500
)




# call peaks using MACS2
peaks <- CallPeaks(pbmc, macs2.path = "/home/toolkit/local/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

#######################################
saveRDS(peaks, './data/peaks_macs2.rds')
####################################

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)


DefaultAssay(pbmc) <- "ATAC"

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments =  Fragments(pbmc),
  annotation = annotations
)


DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)


DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)



allen <- readRDS("./data/allen_brain.rds")

# use the RNA assay in the SNARE-seq data for integration with scRNA-seq
DefaultAssay(pbmc) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = allen,
  query = pbmc,
  dims = 1:30,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen$subclass,
  weight.reduction = pbmc[['pca']],
  dims = 1:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)




DefaultAssay(pbmc) <- "RNA"

pbmc <- FindVariableFeatures(pbmc, nfeatures = 3000)
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 30)
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction.name = "umap")
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 3)



pbmc$type=pbmc$predicted.id


DimPlot(pbmc, group.by = 'type', label = TRUE, reduction = 'umap')+NoLegend()





new.cluster.ids <- c(
  "Neuron.Ex",
  "Neuron.Ex",
  "Neuron.Ex",
  "Neuron.Ex",
  "Neuron.Ex",
  "Neuron.Ex",
  "Neuron.In",
  "Neuron.In",
  "Astro",
  "Neuron.In",
  "ODC",
  "Neuron.Ex",
  "Neuron.Ex",
  "Unclear"
)



Idents(pbmc)=pbmc$seurat_clusters
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids)
pbmc$type <- Idents(pbmc)


p1=DimPlot(pbmc, group.by = 'seurat_clusters', label = TRUE, reduction = 'umap')+NoLegend()
p2=DimPlot(pbmc, group.by = 'predicted.id', label = TRUE, reduction = 'umap')+NoLegend()

p1+p2

DimPlot(pbmc, group.by = 'type', label = TRUE, reduction = 'umap')+NoLegend()


saveRDS(pbmc, file='scATAC_pbmc.rds')






















