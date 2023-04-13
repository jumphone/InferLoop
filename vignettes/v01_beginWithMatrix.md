

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
    
    count_mat = inferloop.loadSignal('count_mat.txt')
    
    
