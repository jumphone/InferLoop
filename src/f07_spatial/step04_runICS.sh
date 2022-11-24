

OLD_NET=./scATAC_loop.txt
NET=./ics_out/UNIQ_LOOP.txt

python ../src00_ics/step0_uniqNet.py $OLD_NET $NET


MAT=./scATAC_normalized_exp_clustered.txt

INDEX=./ics_out/this.index

python3 ../src00_ics/step1_BuildIndex.py $NET $MAT $INDEX 
###########################
TAG=r0
nohup python3 ../src00_ics/step2_CalICS_$TAG\.py   $INDEX ics_out/this.out.$TAG &













