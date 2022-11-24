
#/home/toolkit/tools/R4.2.0/bin/Rscript
#/home/toolkit/tools/R4.2.0/bin/R


source('https://gitee.com/jumphone/VISA/raw/master/VISA.R')

mat=read.table('signal_mat.txt', sep='\t', row.names=1, header=TRUE,check.names=FALSE)
a=t(mat)[,1]
b=t(mat)[,2]

type=readRDS('signal_type.rds')

setwd('/home/disk/database/data/ICS_scHIC/f00_demo')

CEX=2
PCH=15
COL_LABEL=c('blue','skyblue','grey90','red','gold1')
COL_VALUE=c(-2,-0.5,0,0.5,2)


plot(density(a,bw=0.05),type='h',lwd=2,col='royalblue1')


pdf('f00p05_signal_demo.pdf',width=3.5,height=4)
###########################
Z=t(read.table('signal_ics_out/this.out.ori', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
Z1=Z
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('signal_ics_out/this.out.r0', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
Z2=Z
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
#################
dev.off()

pdf('f00p06_signal_boxplot_demo.pdf',width=2.5,height=3.5)
boxplot(Z1[which(type=='OPC')],Z1[which(type!='OPC')],pch='+',col=c('indianred1','grey60'),cex=0.5)
t.test(Z1[which(type=='OPC')],Z1[which(type!='OPC')])
# p-value = 0.007117

boxplot(Z2[which(type=='OPC')],Z2[which(type!='OPC')],pch='+',col=c('indianred1','grey60'),cex=0.5)
t.test(Z2[which(type=='OPC')],Z2[which(type!='OPC')])
#p-value = 0.0006384

dev.off()


