

NET=/home/disk/database/data/ICS_scHIC/f01_frontalCortex/ics_out/UNIQ_LOOP.txt



TAG=N50
MAT=./clst/scATAC_$TAG\_normalized_exp_clustered.txt
INDEX=./clst/scATAC_$TAG\.index
python3 ../src00_ics/step1_BuildIndex.py $NET $MAT $INDEX 
nohup python3 ../src00_ics/step2_CalICS_r0.py $INDEX ./clst/scATAC_$TAG\.out.r0 &


TAG=N100
MAT=./clst/scATAC_$TAG\_normalized_exp_clustered.txt
INDEX=./clst/scATAC_$TAG\.index
python3 ../src00_ics/step1_BuildIndex.py $NET $MAT $INDEX
nohup python3 ../src00_ics/step2_CalICS_r0.py $INDEX ./clst/scATAC_$TAG\.out.r0 &


TAG=N150
MAT=./clst/scATAC_$TAG\_normalized_exp_clustered.txt
INDEX=./clst/scATAC_$TAG\.index
python3 ../src00_ics/step1_BuildIndex.py $NET $MAT $INDEX
nohup python3 ../src00_ics/step2_CalICS_r0.py $INDEX ./clst/scATAC_$TAG\.out.r0 &


TAG=N200
MAT=./clst/scATAC_$TAG\_normalized_exp_clustered.txt
INDEX=./clst/scATAC_$TAG\.index
python3 ../src00_ics/step1_BuildIndex.py $NET $MAT $INDEX
nohup python3 ../src00_ics/step2_CalICS_r0.py $INDEX ./clst/scATAC_$TAG\.out.r0 &















