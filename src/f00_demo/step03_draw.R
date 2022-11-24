
#/home/toolkit/tools/R4.2.0/bin/Rscript
#/home/toolkit/tools/R4.2.0/bin/R


source('https://gitee.com/jumphone/VISA/raw/master/VISA.R')

mat=read.table('mat.txt', sep='\t', row.names=1, header=TRUE,check.names=FALSE)
a=t(mat)[,1]
b=t(mat)[,2]


setwd('/home/disk/database/data/ICS_scHIC/f00_demo')

CEX=1
PCH=15
COL_LABEL=c('blue','skyblue','grey90','red','gold1')
COL_VALUE=c(-2,-0.5,0,0.5,2)


pdf('p01_demo.pdf',width=3.5,height=4)
###########################
Z=t(read.table('ics_out/this.out.mean', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.hmean', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.ori', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.nr200', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.nr150', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.nr100', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.nr50', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.r0', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.r50', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.r100', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.r150', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.r200', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
Z=t(read.table('ics_out/this.out.gmean', sep='\t', row.names=1, header=TRUE,check.names=FALSE))[,1]
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)
###########################
#################
dev.off()


