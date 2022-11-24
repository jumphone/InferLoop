
BEDPE=PDGFRA_TSS_hg19.1k.LOOP.tss.bed.demo.bed
SRC=/home/disk/database/data/scHIC/GSE130711_analysis/src/coolLoopCount.py
OUT=./demo_signal/

TAG=Astro
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool.fithic/matrix_10k.h5.norm.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &
TAG=Neuron.Ex
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool.fithic/matrix_10k.h5.norm.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &
TAG=Neuron.In
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool.fithic/matrix_10k.h5.norm.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &
TAG=MG
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool.fithic/matrix_10k.h5.norm.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &
TAG=ODC
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool.fithic/matrix_10k.h5.norm.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &
TAG=OPC
COOL=/home/disk/database/data/scHIC/GSE130711_cool/Human_cluster_cool/$TAG\_all_brain.txt_1kb_contacts.cool.fithic/matrix_10k.h5.norm.cool
nohup python3 $SRC $COOL $BEDPE $OUT\/$TAG\.coolCount.txt &

