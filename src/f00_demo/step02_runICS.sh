

:<<!
python3 ../src00_ics/step1_BuildIndex.py net.txt mat.txt ics_out/this.index


python3 ../src00_ics/step2_CalICS_mean.py  ics_out/this.index ics_out/this.out.mean
python3 ../src00_ics/step2_CalICS_hmean.py ics_out/this.index ics_out/this.out.hmean
python3 ../src00_ics/step2_CalICS_ori.py   ics_out/this.index ics_out/this.out.ori
!
python3 ../src00_ics/step2_CalICS_gmean.py ics_out/this.index ics_out/this.out.gmean
:<<!
TAG=nr200
python3 ../src00_ics/step2_CalICS_$TAG\.py   ics_out/this.index ics_out/this.out.$TAG
TAG=nr150
python3 ../src00_ics/step2_CalICS_$TAG\.py   ics_out/this.index ics_out/this.out.$TAG
TAG=nr100
python3 ../src00_ics/step2_CalICS_$TAG\.py   ics_out/this.index ics_out/this.out.$TAG
TAG=nr50
python3 ../src00_ics/step2_CalICS_$TAG\.py   ics_out/this.index ics_out/this.out.$TAG
TAG=r0
#python3 ../src00_ics/step2_CalICS_$TAG\.py   ics_out/this.index ics_out/this.out.$TAG
TAG=r50
#python3 ../src00_ics/step2_CalICS_$TAG\.py   ics_out/this.index ics_out/this.out.$TAG
TAG=r100
python3 ../src00_ics/step2_CalICS_$TAG\.py   ics_out/this.index ics_out/this.out.$TAG
TAG=r150
python3 ../src00_ics/step2_CalICS_$TAG\.py   ics_out/this.index ics_out/this.out.$TAG
TAG=r200
python3 ../src00_ics/step2_CalICS_$TAG\.py   ics_out/this.index ics_out/this.out.$TAG
!


