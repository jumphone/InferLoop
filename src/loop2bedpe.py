import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')
for line in fi:
    seq=line.rstrip().split('\t')
    fo.write(seq[0].replace('-','\t')+'\t'+seq[1].replace('-','\t')+'\n')

fo.close()
