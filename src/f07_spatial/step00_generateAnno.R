library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

annotations@seqnames@values=paste0('chr',annotations@seqnames@values)
annotations@seqinfo@seqnames=paste0('chr',annotations@seqinfo@seqnames)
annotations@seqinfo@genome=rep('mm10',length(annotations@seqinfo@genome))
annotations@seqnames@values[which(annotations@seqnames@values=='chrMT')]='chrM'
annotations@seqinfo@seqnames[which(annotations@seqinfo@seqnames=='chrMT')]='chrM'
annotations@seqnames@values=factor(annotations@seqnames@values,levels=annotations@seqinfo@seqnames)

saveRDS(annotations, '/home/database/annotation/mm10/hg10_signac_ucsc_annotations.rds')

