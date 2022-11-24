


.getSub<-function(TAG){

###################
LOOP=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(1))$V1
SSN=fread(paste0("./ics_out/this.out.",TAG),header=F,sep='\t',select=c(2:(ncol(pbmc)+1)))
SSN=as.matrix(SSN)
rownames(SSN)=LOOP
colnames(SSN)=colnames(pbmc)

##########
#CICERO_LOOP=fread(paste0("./ics_out/this.out.cicero"),header=F,sep='\t',select=c(1))$V1
#SSN=SSN[which(rownames(SSN) %in% CICERO_LOOP),]
##########

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


#CCC=apply(COR,2,scale,center=F)
#CCC=COR
CCC=apply(COR,2,rank)
#CCC=COR
rownames(CCC)=rownames(COR)
colnames(CCC)=MATCH[,1]

OUT=c()
i=1
while(i<=ncol(CCC)){
    i_organ=colnames(CCC)[i]
    this_index=which(this_organ %in% i_organ)
    if(length(this_index)>1){
        this_out=apply(CCC[this_index,],2,max)
        }else{
        this_out=CCC[this_index,]
        }
    OUT=cbind(OUT,this_out)
    i=i+1}

colnames(OUT)=colnames(CCC)

TMP=OUT
diag(TMP)=NA

A=diag(OUT)
B=as.numeric(TMP)
B=B[which(!is.na(B))]

#this_sub=mean(A)-mean(B)
return(A)
}





