#  /home/toolkit/tools/R4.2.0/bin/R


library(Signac)
library(Seurat)
library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v86) #hg38
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

opath='/home/disk/database/data/ICS_scHIC/fxx_spatial/'

fpath=paste0('/home/disk/database/data/ICS_scHIC/fxx_spatial/download/GSM5238386_ME13_50um.fragments.tsv.gz')
bpath=paste0('/home/database/reference/mm10/mm10.fa.size.100k_coolBin.bed')


BIN=read.table(bpath,sep='\t',header=F,row.names=NULL)

TMP=BIN
PEAK=GRanges(seqnames = TMP$V1,
            ranges = IRanges(start = TMP$V2+1,
                             end = TMP$V3,
                             names = paste0('bin',c(1:nrow(TMP)))
                             )
             )


fragments <- CreateFragmentObject(fpath,validate.fragments=FALSE)


MAT=FeatureMatrix(
      fragments,
      features=PEAK,
      cells = NULL,
      process_n = 2000,
      sep = c(":", "-"),
      verbose = TRUE
      )

saveRDS(MAT, paste0(opath,'/mat_100k.rds'))

rpath=paste0(opath,'/mat_100k.rds')


MAT=readRDS(rpath)

chrom_assay <- CreateChromatinAssay(
  counts = MAT,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fpath,
  min.cells = 20,
  min.features = 500
  )



pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
  )

################################
library(GenomeInfoDb)
annotations=readRDS('/home/database/annotation/mm10/mm10_signac_ucsc_annotations.rds')

Annotation(pbmc) <- annotations

#pbmc <- RunTFIDF(pbmc)


###########################################
gene.activities <- GeneActivity(pbmc)

pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

#######################################

library(Seurat)
library(Signac)


peaks <- CallPeaks(pbmc, macs2.path = "/home/toolkit/local/bin/macs2")

saveRDS(peaks, file='scATAC_pbmc_peaks_macs2.rds')

peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg19, invert = TRUE)


macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
  )


pbmc[["macs2"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fpath,
  annotation = readRDS('/home/database/annotation/mm10/mm10_signac_ucsc_annotations.rds')
  )



DefaultAssay(pbmc)='peaks'
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)

DimPlot(pbmc)

saveRDS(pbmc, 'scATAC_pbmc.rds')






























































