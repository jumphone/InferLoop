
GENOME=/home/database/reference/hg38/hg38.fa.size
GENE=./eQTL/gencode.v39.annotation.gtf.pc.sorted.changed.bed.uniqGene.sorted.bed
TSS=./eQTL/gencode.v39.annotation.gtf.pc.sorted.changed.bed.uniqGene.sorted.bed.TSS.bed

python ./eQTL/getUniq_TSS.py $TSS $TSS\.uniq.bed 
bedtools sort -i $TSS\.uniq.bed > $TSS\.sorted.bed

TSS=./eQTL/gencode.v39.annotation.gtf.pc.sorted.changed.bed.uniqGene.sorted.bed.TSS.bed.sorted.bed
SNP=./eQTL/ALL_EQTL_1k_merged.bed

bedtools slop -i $GENE -g $GENOME -b 10000 > $GENE\.slop10k.bed
bedtools subtract -a $SNP -b $GENE\.slop10k.bed  > $SNP\.noGENE10k.bed

SNP=./eQTL/ALL_EQTL_1k_merged.bed.noGENE10k.bed
bedtools window -a $SNP -b $TSS -w 500000 > ./eQTL/SNP_TSS_500k.bed


SNP_TSS=./eQTL/SNP_TSS_500k.bed
SLOP=1000
python bedpe2regionLoop.py $SNP_TSS $SNP_TSS\.region.bed $SNP_TSS\.loop.txt $SLOP


