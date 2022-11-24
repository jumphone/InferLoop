import sys
ALL_LOOP=sys.argv[1]
L1SNP=sys.argv[2]
L1TSS=sys.argv[3]
L2SNP=sys.argv[4]
L2TSS=sys.argv[5]
fo=open(sys.argv[6],'w')

LOOP={}
fi=open(ALL_LOOP)
i=1
for line in fi:
    seq=line.rstrip().split('\t')
    LOOP['L'+str(i)]='\t'.join(seq[0:6])
    i=i+1


f1=open(L1SNP)
f2=open(L2TSS)
SNP={}
TSS={}
for line in f1:
    seq=line.rstrip().split('\t')
    this_lp=seq[3].split('P')[0]
    #print(this_lp)
    this_snp='\t'.join(seq[4:])
    SNP[this_lp]=this_snp

for line in f2:
    seq=line.rstrip().split('\t')
    this_lp=seq[3].split('P')[0]
    this_tss='\t'.join(seq[4:])
    TSS[this_lp]=this_tss

for one in SNP:
    if one in TSS:
        fo.write(LOOP[one]+'\t'+one +'\t'+TSS[one]+'\t'+SNP[one]+'\t'+'snp_tss'+'\n')


f1=open(L2SNP)
f2=open(L1TSS)
SNP={}
TSS={}
for line in f1:
    seq=line.rstrip().split('\t')
    this_lp=seq[3].split('E')[0]
    this_snp='\t'.join(seq[4:])
    SNP[this_lp]=this_snp

for line in f2:
    seq=line.rstrip().split('\t')
    this_lp=seq[3].split('E')[0]
    this_tss='\t'.join(seq[4:])
    TSS[this_lp]=this_tss

for one in SNP:
    if one in TSS:
        fo.write(LOOP[one]+'\t'+one +'\t'+TSS[one]+'\t'+SNP[one]+'\t'+'tss_snp'+'\n')





