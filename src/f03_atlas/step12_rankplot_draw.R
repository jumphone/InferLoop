
#/home/toolkit/tools/R4.2.0/bin/R

library(Seurat)
library(Signac)
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
library(data.table)
library(stringr)



pbmc=readRDS('./eQTL/scAtlas_pbmc.rds')
SNP_TSS=read.table('./eQTL/SNP_TSS_500k.bed.loop.txt',sep='\t',row.names=NULL,header=FALSE)
SNP_TSS_TAG=SNP_TSS[,3]
SNP_TSS_GENE=SNP_TSS[,4]
INFO=cbind(paste0(SNP_TSS[,1],'.And.',SNP_TSS[,2]),SNP_TSS_TAG,SNP_TSS_GENE)
rownames(INFO)=INFO[,1]

PAIR=paste0(INFO[,2],'_',INFO[,3])

MATCH=read.table('MATCH.txt',sep='\t')


#############################



TAG='r0'


###################
LOOP=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(1))$V1
SSN=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(2:(ncol(pbmc)+1)))
SSN=as.matrix(SSN)
rownames(SSN)=LOOP
colnames(SSN)=colnames(pbmc)
###########
dim(SSN)
dim(INFO)

D1=SSN[which(rownames(SSN) %in% rownames(INFO)),]
D2=INFO[which(rownames(INFO) %in% rownames(SSN)),]
D1=D1[order(rownames(D1)),]
D2=D2[order(rownames(D2)),]
#set.seed(123)
#DD1=apply(D1,2,function(x){ecdf(x)(x)})


this_tag=D2[,2]
this_gene=D2[,3]
#this_organ=colnames(D1)
this_organ=t(matrix(unlist(str_split(colnames(D1),'_')),nrow=2))[,1]
this_pair=paste0(this_tag,'_',this_gene)



COR=c()
i=1
while(i<=nrow(MATCH)){

print(i)
    eq_file=paste0('./eQTL/GTEx_Analysis_v8_eQTL/',MATCH[i,2],'/tag.bed')
    eq_data=as.matrix(fread(eq_file))
    eq_tag=paste0(eq_data[,1],'x',eq_data[,2],'x',eq_data[,3])
    eq_tag=str_replace_all(eq_tag,' ','')
    eq_gene=eq_data[,4]
    eq_pair=paste0(eq_tag,'_',eq_gene)
    #eq_lpv=-log10(as.numeric(eq_data[,5]))
    eq_lpv=as.numeric(eq_data[,5])


X1=D1
rownames(X1)=this_pair
X2=cbind(eq_lpv,eq_lpv)
rownames(X2)=eq_pair
X1=X1[which(rownames(X1) %in% rownames(X2)),]
X2=X2[which(rownames(X2) %in% rownames(X1)),]
X1=X1[order(rownames(X1)),]
X2=X2[order(rownames(X2)),]


this_cor=cor(X1,X2[,1])
COR=cbind(COR,this_cor)
i=i+1
}



colnames(COR)=MATCH[,1]


CCC=apply(COR,2,rank)
#CCC=COR
rownames(CCC)=rownames(COR)
colnames(CCC)=MATCH[,1]


write.table(CCC,file='eQTL_tissue_type_rank_table.txt',sep='\t',row.names=T,col.names=T,quote=F)



OUT=CCC
OUT=OUT*0
TMP=c()
i=1
while(i<=ncol(CCC)){
    i_organ=colnames(CCC)[i]
    this_index=which(this_organ %in% i_organ)
    ###################
    this_out=max(CCC[this_index,i])
    print(this_out)
    ###########################
    TMP=c(TMP,this_out)
    OUT[this_index,i]=1
    OUT[which(CCC[,i]==this_out),i]=2
    OUT[,i]=OUT[order(-CCC[,i]),i]
    i=i+1}


OUT=OUT[,order(-TMP)]
rownames(OUT)=c(nrow(OUT):1)




library('ComplexHeatmap')
library('circlize')

#mat=mat
o.mat=OUT
col_fun =colorRamp2(c(0,1,2 ), c('white','grey50','red'))
ht=Heatmap(o.mat,row_title='',name="rank",cluster_rows=F,cluster_columns=F,
        show_column_dend = F, show_row_dend = F,
        show_column_names=FALSE, show_row_names=FALSE,
        col=col_fun, border = TRUE,column_split=factor(colnames(o.mat),levels=colnames(o.mat)),
         column_title=' '
        )

ht


pdf('p13_ICS_r0_heatmap.pdf',width=4,height=7)
print(ht)
dev.off()







































