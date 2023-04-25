    ...

    TMP=getMatrixFromProject(proj,useMatrix ='PeakMatrix')
    CHR=TMP@rowRanges@seqnames
    START=TMP@rowRanges@ranges@start
    END=TMP@rowRanges@ranges@start+TMP@rowRanges@ranges@width-1
    PEAK_MAT=TMP@assays@data$PeakMatrix
    colnames(PEAK_MAT)=stringr::str_replace_all(rownames(TMP@colData),'#','_')
    rownames(PEAK_MAT)=paste0(CHR,'-',START,'-',END)

    library(Seurat)
    library(Signac)
    peak_mat=CreateSeuratObject(counts=PEAK_MAT,assay='peak')
    peak_mat=RunTFIDF(peak_mat)
    
    source('InferLoop.R')
    coa=getCoAccessibility(
        ArchRProj = proj,
        corCutOff = 0.1,
        resolution = 1,
        returnLoops = FALSE
        )

    bed=cbind(
              as.character(coa@metadata$peakSet@seqnames),
              coa@metadata$peakSet@ranges@start,
              coa@metadata$peakSet@ranges@start+coa@metadata$peakSet@ranges@width-1
             )

    bedpe=cbind(bed[coa$queryHits,],bed[coa$subjectHits,])
    value=coa$correlation

    conns=data.frame(
               p1=paste0(bedpe[,1],'-',bedpe[,2],'-',bedpe[,3]),
               p2=paste0(bedpe[,4],'-',bedpe[,5],'-',bedpe[,6]),
               cor=value
               )

    inferloop.writeNet(conns, "net.txt", cut=400000)
    write.table(peak_mat[['peak']]@data,file='mat.txt', row.names=T,col.names=T,quote=F,sep='\t')


    source('InferLoop.R')
    mat=inferloop.loadSignal('mat.txt')
    net=read.table('net.txt',sep='\t',header=F,row.names=NULL)
    net_uniq=inferloop.getUniqLoop(net)
    MAT=inferloop.inferLoopSignal(mat, net_uniq, r=0)
    write.table(MAT,file='signal_mat.txt', row.names=T,col.names=T,quote=F,sep='\t')
