
TAG='400k'



# /home/toolkit/tools/R4.2.0/bin/R

library(Signac)
library(GenomicRanges)
library(Seurat)
#args <- commandArgs(TRUE)
#opath=args[1]
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')


COL_LIST=list(DType=c('scATAC'='indianred1','scHIC'='royalblue1'),
              CType=c('Astro'='red','MG'='gold1','Neuron.Ex'='darkgreen','Neuron.In'='darkolivegreen1','ODC'='blue3','OPC'='skyblue'))



.draw_ht<-function(mat=mat){

    library('ComplexHeatmap')
    library('circlize')

    DType=t(matrix(unlist(stringr::str_split(colnames(mat),'_')),nrow=2))[,2]
    CType=t(matrix(unlist(stringr::str_split(colnames(mat),'_')),nrow=2))[,1]
    ha_top = HeatmapAnnotation( DType=DType, CType=CType,
                            col=COL_LIST
                          )
    ha_left = rowAnnotation( DType=DType, CType=CType,
                         col=COL_LIST
                          )
    col_fun =colorRamp2(c(-1,-0.2,-0.1,0,0.05,0.1,0.2,1 ), c('blue','skyblue','grey80','grey95','grey80','indianred1','red1','gold1'))

    mat=mat
    o.mat=t(mat)
    ht=Heatmap(o.mat,row_title='',name="cor",cluster_rows=TRUE,cluster_columns=TRUE,
        show_column_dend = TRUE, show_row_dend = TRUE,
        show_column_names=FALSE, show_row_names=FALSE,
        col=col_fun, border = TRUE,
        row_names_side='right',
        top_annotation = ha_top,
        left_annotation = ha_left,
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        #cell_fun = function(j, i, x, y, width, height, fill) {
        #  grid.text(sprintf("%.2f", o.mat[i, j]), x, y, gp = gpar(fontsize = 10))
        #  }
        )
    return(ht)
    }


real.ref=readRDS(paste0('./coolCount_',TAG,'/COM_real.ref.rds'))
real.ref=real.ref[which(rowSums(real.ref)>0),]

#pdf('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/p01_cor_heatmap_all_small.pdf',width=4,height=3.5)
pdf(paste0('./f10p01_cor_heatmap_all_cortex_small_',TAG,'.pdf'),width=5,height=4)
################################################################
#################
pred.ref=readRDS(paste0('./ics_out/COM_pred.ref_',TAG,'.rds'))
COM=.simple_combine(real.ref,pred.ref)
D1=COM$exp_sc_mat1
D2=COM$exp_sc_mat2
ALL=cbind(D1,D2)
COR=cor(ALL,method='pearson')
ht=.draw_ht(COR)
print(ht)
#################


dev.off()















