TSS=/home/database/annotation/hg19/Homo_sapiens.GRCh37.75.chr.pc.gene.SYM.bed.TSS.bed
LOOP=/home/disk/database/data/scHIC/GSE130711_analysis/scATAC_RefLoop.bed
GENOME=/home/database/reference/hg19/hg19.fa.size

grep PDGFRA $TSS > PDGFRA_TSS_hg19.bed

bedtools slop -i PDGFRA_TSS_hg19.bed -g $GENOME -b 1000 > PDGFRA_TSS_hg19.1k.bed 

bedtools pairtobed -a $LOOP -b PDGFRA_TSS_hg19.1k.bed > PDGFRA_TSS_hg19.1k.LOOP.bed 


