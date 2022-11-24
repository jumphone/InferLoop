cat /home/disk/database/data/ICS_scHIC/f03_atlas/eQTL/GTEx_Analysis_v8_eQTL/*format/pair.bed | cut -f 1,2,3 |  bedtools sort -i - | uniq - >  ALL_EQTL.bed

GENOME=/home/database/reference/hg38/hg38.fa.size
SLOP=1000
bedtools slop -i ALL_EQTL.bed -g $GENOME -b $SLOP | bedtools merge -i - > ALL_EQTL_1k_merged.bed

bedtools intersect -a ALL_EQTL_1k_merged.bed -b ALL_EQTL.bed -wa -wb > ALL_EQTL_1k_merged.bed.detail.bed


