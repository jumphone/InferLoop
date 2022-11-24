import sys

fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.TSS.bed','w')


for line in fi:
    seq=line.rstrip().split('\t')
    if seq[5]=='+':
       loc=int(seq[1])+1
    if seq[5]=='-':
       loc=int(seq[2])
    this_out=[seq[0],str(loc-1),str(loc),seq[3],seq[4],seq[5]]
    fo.write('\t'.join(this_out)+'\n')

fi.close()
fo.close()


