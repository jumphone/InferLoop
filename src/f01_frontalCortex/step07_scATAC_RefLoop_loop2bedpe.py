

fi=open('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/tracks/scATAC_RefLoop.txt')


header=fi.readline().rstrip().split('\t')

import os
os.chdir('/home/disk/database/data/ICS_scHIC/f01_frontalCortex/tracks')

fo_list=[]
for one in header:
    fo_list.append(open(one+'.bedpe','w'))



for line in fi:
    seq=line.rstrip().split('\t')
    if 1==1:
        sss=seq[0].split('.And.')
        out=sss[0].replace('-','\t')
        tmp=sss[1].split('-')
        out=out+'\t'+tmp[0]+'\t'+tmp[1]+'\t'+tmp[2]+'\t'
        i=0
        while i <len(fo_list): 
            this_score=round(float(seq[1+i]),3)
            fo=fo_list[i]
            if this_score > 0:
                fo.write(out+str(this_score*100)+'\n')
            i=i+1

for fo in fo_list:
    fo.close()


import subprocess as sp
for one in header:
    sp.Popen('bedtools sort -i '+one+'.bedpe'+' > '+one+'.sorted.bedpe',shell=True).wait()
    sp.Popen('rm -rf '+one+'.bedpe',shell=True).wait()

fi.close()
    


