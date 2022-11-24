
GENOME=/home/database/reference/hg38/hg38.fa.size

TSS=/home/database/annotation/hg38/gencode.v39.annotation.gtf.pc.sorted.bed.TSS.bed
LOOP=/home/disk/database/data/scHIC/Blood/analysis/scATAC_RefLoop.bed


SLOP=1000
bedtools slop -i $TSS -g $GENOME -b $SLOP > $TSS\.slop$SLOP\.bed


bedtools pairtobed -a $LOOP -b $TSS\.slop$SLOP\.bed > RNA/LoopAndTSS.bed 


