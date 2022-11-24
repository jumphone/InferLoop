

NET=./scATAC_loop.txt
python3 ../src00_ics/step0_uniqNet.py $NET $NET\.uniq.txt


NET=./scATAC_loop.txt.uniq.txt
MAT=./scATAC_normalized_exp_clustered.txt
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



