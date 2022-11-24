
import subprocess
import multiprocessing
import os


RS='/home/disk/database/data/ICS_scHIC/f03_atlas/single/rs02_getAgg.R'

def Work(input_path,a):
    print input_path
    subprocess.Popen('/home/toolkit/tools/R4.2.0/bin/Rscript '+RS+' '+input_path,shell=True).wait()


fa=open('/home/disk/database/data/ICS_scHIC/f03_atlas/LIST.txt')
PROC_LIMIT=30
jobs=[]
i=1

for line in fa:
    print i;i+=1
    seq=line.rstrip().split('\t')
    input_path='/home/disk/database/data/scATAC_ATLAS/raw/FRAG/'+seq[0]+'.format'
    print input_path
    if 1==1:
        p=multiprocessing.Process(target=Work, args=(input_path,1))
        p.start()
        jobs.append(p)
        if len(jobs)>=PROC_LIMIT:
            for p in jobs:
                p.join()
            jobs=[]
for p in jobs:
    p.join()



