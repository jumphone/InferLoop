
GENOME=/home/database/reference/mm10/mm10.fa.size

TSS=/home/database/annotation/mm10/Mus_musculus.GRCm38.87.chr.gtf.pc.bed.TSS.sorted.bed
LOOP=./scATAC_RefLoop.bed


SLOP=1000
bedtools slop -i $TSS -g $GENOME -b $SLOP > $TSS\.slop$SLOP\.bed


bedtools pairtobed -a $LOOP -b $TSS\.slop$SLOP\.bed > RNA/LoopAndTSS.bed 


