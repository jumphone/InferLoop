

python3 ../src00_ics/step1_BuildIndex.py signal_net.txt signal_mat.txt signal_ics_out/this.index


python3 ../src00_ics/step2_CalICS_ori.py   signal_ics_out/this.index signal_ics_out/this.out.ori

TAG=r0
python3 ../src00_ics/step2_CalICS_$TAG\.py   signal_ics_out/this.index signal_ics_out/this.out.$TAG


TAG=r50
python3 ../src00_ics/step2_CalICS_$TAG\.py   signal_ics_out/this.index signal_ics_out/this.out.$TAG

TAG=r100
python3 ../src00_ics/step2_CalICS_$TAG\.py   signal_ics_out/this.index signal_ics_out/this.out.$TAG
