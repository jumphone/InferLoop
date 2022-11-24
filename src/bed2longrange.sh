BED=$1

python3 /home/database/data/HiChip_Li/EACH/bed2longrange.py $BED 
bedtools sort -i $BED\.longrange > $BED\.longrange.bed
bgzip $BED\.longrange.bed
tabix -p bed $BED\.longrange.bed.gz





