
import subprocess
import multiprocessing
import os


SCRIPT='/home/disk/database/data/ICS_scHIC/f03_atlas/eQTL/src/step05_getUniq.py'

def Work(input_path,a):
    print input_path
    subprocess.Popen('python '+SCRIPT+' '+input_path+'.format/tag.bed '+input_path+'.format/tag_uniq.bed',shell=True).wait()


fa=open('/home/disk/database/data/ICS_scHIC/f03_atlas/eQTL/LIST.txt')
PROC_LIMIT=1
jobs=[]
i=1

for line in fa:
    print i;i+=1
    seq=line.rstrip().split('\t')
    input_path='/home/disk/database/data/ICS_scHIC/f03_atlas/eQTL/GTEx_Analysis_v8_eQTL/'+seq[0]
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



