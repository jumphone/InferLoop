:<<!
tabix -p bed GSM6043255_ME11_20um_fragments.tsv.gz 

samtools faidx /home/database/reference/mm10/mm10.fa

cut -f 1,2 /home/database/reference/mm10/mm10.fa.fai > /home/database/reference/mm10/mm10.fa.size
!

CHROM_SIZE=/home/database/reference/mm10/mm10.fa.size

BINSIZE=10000
cooler makebins $CHROM_SIZE $BINSIZE > $CHROM_SIZE\.10k_coolBin.bed

BINSIZE=100000
cooler makebins $CHROM_SIZE $BINSIZE > $CHROM_SIZE\.100k_coolBin.bed
