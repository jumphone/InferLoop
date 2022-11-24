
a=read.table('PDGFRA_TSS_hg19.1k.LOOP.tss.bed.demo.bed.signal',sep='\t')

b=c(as.numeric(a[2,12]),mean(as.numeric(a[2,7:11])))

pdf('f00p04_demo_pdgfra.pdf',width=5,height=3)
barplot(as.numeric(a[2,7:12]),col=c(rep('grey60',5),'red'),border=F)
abline(h=1)
dev.off()


