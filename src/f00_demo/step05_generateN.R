
i=1
while(i<=10){
N=i*1000
set.seed(N)
a=runif(N)
b=runif(N)

out=rbind(a,b)
rownames(out)=c('a','b')
colnames(out)=c(1:ncol(out))
write.table(out,  file=paste0('./testN/mat',N,'.txt'), row.names=T,col.names=T,quote=F,sep='\t')

network=matrix(NA, nrow=1,ncol=2)
network[1,1]='a'
network[1,2]='b'
write.table(network,  file=paste0('./testN/net',N,'.txt'), row.names=F,col.names=F,quote=F,sep='\t')

i=i+1
}


