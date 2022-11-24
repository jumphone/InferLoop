

NET=./eQTL/SNP_TSS_500k.bed.loop.txt
MAT=./eQTL/scAtlas_normalizedData.txt

INDEX=./ics_out/this.index

python3 ../src00_ics/step1_BuildIndex.py $NET $MAT $INDEX 


nohup python3 ../src00_ics/step2_CalICS_mean.py  $INDEX ics_out/this.out.mean &
nohup python3 ../src00_ics/step2_CalICS_hmean.py  $INDEX ics_out/this.out.hmean &
nohup python3 ../src00_ics/step2_CalICS_ori.py  $INDEX ics_out/this.out.ori &
###########################
TAG=r0
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &
TAG=r50
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &
#


