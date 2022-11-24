

#/home/toolkit/tools/R4.2.0/bin/Rscript
#/home/toolkit/tools/R4.2.0/bin/R


source('/home/disk/database/data/ICS_scHIC/src00_ics/ICS.R')
source('https://gitee.com/jumphone/VISA/raw/master/VISA.R')


set.seed(123)

N=50

a=rep(1:N,each=N)/N
b=rep(1:N,times=N)/N


N=1000
a=rnorm(N)
b=a*0.3 + rnorm(N)*0.5



CEX=1
PCH=15
COL_LABEL=c('blue','skyblue','grey90','red','gold1')
COL_VALUE=c(-2,-0.5,0,0.5,2)
Z=calICS(a,b,0,4)
###########################
COL=visa.vcol(Z,COL_VALUE,COL_LABEL)
plot(a,b,col=COL,pch=PCH,cex=CEX)




