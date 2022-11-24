BED=$1

python3 bed2longrange.py $BED 
bedtools sort -i $BED\.longrange > $BED\.longrange.bed
bgzip $BED\.longrange.bed
tabix -p bed $BED\.longrange.bed.gz





