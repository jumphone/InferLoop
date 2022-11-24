
set.seed(123)

N=50

a=rep(1:N,each=N)/N
b=rep(1:N,times=N)/N



out=rbind(a,b)
rownames(out)=c('a','b')
colnames(out)=c(1:ncol(out))
write.table(out,  file='mat.txt', row.names=T,col.names=T,quote=F,sep='\t')

network=matrix(NA, nrow=1,ncol=2)
network[1,1]='a'
network[1,2]='b'
write.table(network,  file='net.txt', row.names=F,col.names=F,quote=F,sep='\t')




