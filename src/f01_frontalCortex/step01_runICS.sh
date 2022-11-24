

OLD_NET=/home/disk/database/data/scHIC/GSE130711_analysis/scATAC_loop.txt
NET=./ics_out/UNIQ_LOOP.txt

python ../src00_ics/step0_uniqNet.py $OLD_NET $NET


MAT=/home/disk/database/data/scHIC/GSE130711_analysis/scATAC_normalized_exp_clustered.txt

INDEX=/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/this.index

python3 ../src00_ics/step1_BuildIndex.py $NET $MAT $INDEX 

nohup python3 ../src00_ics/step2_CalICS_mean.py  $INDEX ics_out/this.out.mean &
nohup python3 ../src00_ics/step2_CalICS_hmean.py  $INDEX ics_out/this.out.hmean &
nohup python3 ../src00_ics/step2_CalICS_ori.py  $INDEX ics_out/this.out.ori &

###########################

TAG=nr200
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &
TAG=nr150
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &
TAG=nr100
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &
TAG=nr50
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &
TAG=r0
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &
TAG=r50
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &
TAG=r100
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &
TAG=r150
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &
TAG=r200
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &













