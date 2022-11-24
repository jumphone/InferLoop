
MAT=/home/disk/database/data/scHIC/GSE130711_analysis/scATAC_normalized_exp_clustered.txt

TAG=100k
NET=./scATAC_loop_$TAG\.txt
INDEX=./ics_out/scATAC_$TAG\.index
#python3 ../src00_ics/step0_uniqNet.py $NET $NET\.uniq 
nohup python3 loop2bedpe.py $NET\.uniq $NET\.bed &
#python3 ../src00_ics/step1_BuildIndex.py $NET\.uniq $MAT $INDEX 
#nohup python3 ../src00_ics/step2_CalICS_r0.py $INDEX ./ics_out/scATAC_$TAG\.out.r0 &

TAG=200k
NET=./scATAC_loop_$TAG\.txt
INDEX=./ics_out/scATAC_$TAG\.index
#python3 ../src00_ics/step0_uniqNet.py $NET $NET\.uniq 
nohup python3 loop2bedpe.py $NET\.uniq $NET\.bed &
#python3 ../src00_ics/step1_BuildIndex.py $NET\.uniq $MAT $INDEX 
#nohup python3 ../src00_ics/step2_CalICS_r0.py $INDEX ./ics_out/scATAC_$TAG\.out.r0 &

TAG=300k
NET=./scATAC_loop_$TAG\.txt
INDEX=./ics_out/scATAC_$TAG\.index
#python3 ../src00_ics/step0_uniqNet.py $NET $NET\.uniq 
nohup python3 loop2bedpe.py $NET\.uniq $NET\.bed &
#python3 ../src00_ics/step1_BuildIndex.py $NET\.uniq $MAT $INDEX 
#nohup python3 ../src00_ics/step2_CalICS_r0.py $INDEX ./ics_out/scATAC_$TAG\.out.r0 &

TAG=400k
NET=./scATAC_loop_$TAG\.txt
INDEX=./ics_out/scATAC_$TAG\.index
#python3 ../src00_ics/step0_uniqNet.py $NET $NET\.uniq 
nohup python3 loop2bedpe.py $NET\.uniq $NET\.bed &
#python3 ../src00_ics/step1_BuildIndex.py $NET\.uniq $MAT $INDEX 
#nohup python3 ../src00_ics/step2_CalICS_r0.py $INDEX ./ics_out/scATAC_$TAG\.out.r0 &


















