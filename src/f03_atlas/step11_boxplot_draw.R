

#/home/toolkit/tools/R4.2.0/bin/R



#A=readRDS('compare_real.rds')
A=readRDS('compare_real_scaleCicero.rds')

B=readRDS('compare_random.rds')


a=apply(A,2,median)
b=apply(B,2,median)


length(which(b<a[4]))/500

pdf('p11_compare_boxplot.pdf',width=4,height=5)
boxplot(A,col=c(rep('grey60',4),rep('indianred1',2)),pch='+',ylim=c(0,160),lwd=1,
at = seq(1, ncol(A) * 2, by = 2)  )
abline(h=160)
dev.off()

db=density(b,bw=0.5)

pdf('p12_compare_density.pdf',width=5,height=3)
plot(db$x,db$y,type='h',col='grey60',lwd=2,xlim=c(0,160))
points(db$x,db$y,col='grey10',lwd=2,type='l')
abline(v=148.5)
abline(v=160)
dev.off()


pdf('p14_compare_boxplot_onlyICSr0.pdf',width=1.5,height=5)
boxplot(A[,6],col=c(rep('indianred1',1)),pch='+',ylim=c(0,160),lwd=1)
abline(h=160)
dev.off()




