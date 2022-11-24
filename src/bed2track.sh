BED=$1
cut -f 1,2,3 $BED > $BED\.tmp.bed
bgzip -c $BED\.tmp.bed  > $BED\.gz
tabix -p bed $BED\.gz
