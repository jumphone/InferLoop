
COR=c()
i=1
while(i<=10){
N=i*1000
A=read.table(paste0('./testN/this.out.ori.',N),header=T,row.names=1,sep='\t')
B=read.table(paste0('./testN/this.out.r100.',N),header=T,row.names=1,sep='\t')

A=t(A)[,1]
B=t(B)[,1]
this_cor=cor(A,B)
COR=c(COR,this_cor)
i=i+1
}


TIME=read.table('TIME.txt',sep='\t',header=F,row.names=NULL)
TIME=t(TIME)

pdf('p02_time.pdf')
plot(TIME[,1],type='b',lwd=3,pch=16)
points(TIME[,2],type='b',lwd=3,pch='+',col='red')
dev.off()

