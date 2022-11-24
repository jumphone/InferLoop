
import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')

for line in fi:
    seq=line.rstrip().split('\t')
    seq[4]=str(int(seq[4])-1000)
    seq[5]=str(int(seq[5])+1000) 
    fo.write('\t'.join(seq)+'\n')


fo.close()
fi.close()


