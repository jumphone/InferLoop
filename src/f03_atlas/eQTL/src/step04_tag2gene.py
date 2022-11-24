fa=open('/home/disk/database/data/ICS_scHIC/f03_atlas/eQTL/ALL_EQTL_1k_merged.bed.detail.bed')

snp2tag={}

for line in fa:
    seq=line.rstrip().split('\t')
    this_tag=seq[0]+'\t'+seq[1]+'\t'+seq[2]
    this_snp=seq[3]+'\t'+seq[4]+'\t'+seq[5]
    snp2tag[this_snp]=this_tag

import numpy as np

from scipy.stats import combine_pvalues
import sys
fi=open(sys.argv[1])
PV={}
for line in fi:
    seq=line.rstrip().split('\t')
    this_snp=seq[0]+'\t'+seq[1]+'\t'+seq[2]
    this_tag=snp2tag[this_snp]
    this_gene=seq[3]
    #this_score=-np.log10(float(seq[6]))
    this_score=float(seq[6])
    this_out=this_tag+'\t'+this_gene
    if this_out not in PV:
        PV[this_out] = [this_score]
    else:
        PV[this_out].append(this_score)



fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')
old=set()
for line in fi:
    seq=line.rstrip().split('\t')
    this_snp=seq[0]+'\t'+seq[1]+'\t'+seq[2]
    this_tag=snp2tag[this_snp]
    this_gene=seq[3]
    this_out=this_tag+'\t'+this_gene
    if this_out not in old: 
        old.add(this_out)
        #this_cp=combine_pvalues(PV[this_out],method='mudholkar_george')[1]
        #fo.write(this_out+'\t'+str(-np.log10(this_cp))+'\n')
        #print(PV[this_out])
        fo.write(this_out+'\t'+str(np.median(-np.log10(PV[this_out])))+'\n')

fo.close()



