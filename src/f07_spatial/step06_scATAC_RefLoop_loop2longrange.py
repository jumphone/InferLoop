

fi=open('scATAC_RefLoop.txt')


header=fi.readline().rstrip().split('\t')

import os
os.chdir('./scATAC_RefLoop')

fo_list=[]
for one in header:
    fo_list.append(open(one+'.longrange','w'))



for line in fi:
    seq=line.rstrip().split('\t')
    if 1==1:
        sss=seq[0].split('.And.')
        out=sss[0].replace('-','\t')
        tmp=sss[1].split('-')
        out=out+'\t'+tmp[0]+':'+tmp[1]+'-'+tmp[2]+','
        i=0
        while i <len(fo_list): 
            this_score=round(float(seq[1+i]),3)
            fo=fo_list[i]
            if this_score > 0:
                fo.write(out+str(this_score*100)+'\n')
            i=i+1

for fo in fo_list:
    fo.close()



#bedtools sort -i $BED\.longrange > $BED\.longrange.bed
#bgzip $BED\.longrange.bed
#tabix -p bed $BED\.longrange.bed.gz



import subprocess as sp
for one in header:
    sp.Popen('bedtools sort -i '+one+'.longrange'+' > '+one+'.longrange.bed',shell=True).wait()
    sp.Popen('bgzip '+one+'.longrange.bed',shell=True).wait()
    sp.Popen('tabix -p bed '+one+'.longrange.bed.gz',shell=True).wait()
    sp.Popen('rm -rf '+one+'.longrange',shell=True).wait()


fi.close()
    


