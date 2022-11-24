
#TAG=100k
#TAG=200k
#TAG=300k
TAG=400k

SRC=/home/disk/database/data/scHIC/GSE130711_analysis/src/coolLoopCount.py
BEDPE=./scATAC_loop_$TAG\.txt.bed
OUT=./coolCount_$TAG\/
mkdir $OUT
TAG=Astro
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &
TAG=Neuron.Ex
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &
TAG=Neuron.In
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &
TAG=MG
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &
TAG=ODC
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &
TAG=OPC
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &


