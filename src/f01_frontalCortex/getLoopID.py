import sys

fi=open(sys.argv[1])
fo1=open(sys.argv[2]+'.1.bed','w')
fo2=open(sys.argv[2]+'.2.bed','w')

i=1
for line in fi:
    seq=line.rstrip().split('\t')
    out1=[seq[0],seq[1],seq[2],'L'+str(i)+'P1']
    out2=[seq[3],seq[4],seq[5],'L'+str(i)+'P2']
    fo1.write('\t'.join(out1)+'\n')
    fo2.write('\t'.join(out2)+'\n')
    i=i+1

fo1.close()
fo2.close()

