

N=10000
python3 ../src00_ics/step1_BuildIndex.py ./testN/net$N\.txt ./testN/mat$N\.txt ./testN/this.index.$N

echo $[$(date +%s%N)/1000000]
TAG=r100
python3 ../src00_ics/step2_CalICS_$TAG\.py   ./testN/this.index.$N ./testN/this.out.$TAG\.$N
echo $[$(date +%s%N)/1000000]
TAG=ori
python3 ../src00_ics/step2_CalICS_$TAG\.py   ./testN/this.index.$N ./testN/this.out.$TAG\.$N
echo $[$(date +%s%N)/1000000]
