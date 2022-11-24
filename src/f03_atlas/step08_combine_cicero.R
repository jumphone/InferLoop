


#/home/toolkit/tools/R4.2.0/bin/R



library(Seurat)
library(stringr)

LIST=read.table('LIST.txt',sep='\t')

pbmc=readRDS('./eQTL/scAtlas_pbmc.rds')
#ORGAN=t(matrix(unlist(str_split(colnames(pbmc),'_')),nrow=2))[,1]
#TYPE=colnames(pbmc)#t(matrix(unlist(str_split(colnames(pbmc),'_')),nrow=2))[,2]
LOOP=read.table('./eQTL/SNP_TSS_500k.bed.loop.txt')

ORGAN=LIST[,2]


OUT=c()

i=1
while(i<=length(ORGAN)){
    this_organ=ORGAN[i]
    this_file=paste0('/home/disk/database/data/scATAC_ATLAS/raw/FRAG/',LIST[which(LIST[,2]==this_organ),1],'.format/cicero_celltype_loop.rds')
    this_cicero=readRDS(this_file)
    i_out=c()
    j=1
    while(j<=length(this_cicero)){
        this_out=this_cicero[[j]][,3]
        i_out=cbind(i_out, this_out)
        j=j+1
        }
    OUT=cbind(OUT, i_out)
    colnames(OUT)[(ncol(OUT)-length(this_cicero)+1):ncol(OUT)]=names(this_cicero)
    print(i)
    i=i+1
    }


this_tag=paste0(str_replace_all(this_cicero[[1]][,1],'_','-'),'.And.',str_replace_all(this_cicero[[1]][,2],'_','-'))
rownames(OUT)=this_tag


NEW_OUT=OUT
colnames(NEW_OUT)=colnames(pbmc)
i=1
while(i<=ncol(NEW_OUT)){
    this_new_name=colnames(pbmc)[i]
    NEW_OUT[,i]=OUT[,which(colnames(OUT)==this_new_name)]
    print(i)
    i=i+1}

ALL_LOOP=paste0(LOOP[,1],'.And.',LOOP[,2])
USED_OUT=NEW_OUT[which(rownames(NEW_OUT) %in% ALL_LOOP),]
#OTHER=ALL_LOOP[which(!ALL_LOOP %in% rownames(USED_OUT))]
#ADD_OUT=matrix(0,nrow=length(OTHER),ncol=ncol(USED_OUT))
#rownames(ADD_OUT)=OTHER
#colnames(ADD_OUT)=colnames(USED_OUT)
#FINAL_OUT=rbind(USED_OUT,ADD_OUT)



saveRDS(NEW_OUT, file='./ics_out/this.out.cicero_all.rds')


write.table(USED_OUT, file='./ics_out/this.out.cicero',sep='\t',quote=F,row.names=T,col.names=T)




tmp=read.table('./ics_out/this.out.cicero',sep='\t',row.names=1,header=T)

stmp=t(apply(tmp,1,scale))
rownames(stmp)=rownames(tmp)
colnames(stmp)=colnames(tmp)
stmp[which(is.na(stmp))]=0

write.table(stmp, file='./ics_out/this.out.cicero_scale',sep='\t',quote=F,row.names=T,col.names=T)









