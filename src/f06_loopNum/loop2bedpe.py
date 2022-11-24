

import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')

#fi.readline()
i=1
for line in fi:
    seq=line.rstrip().split('\t')
    p1=seq[0]
    p2=seq[1]
    out=p1.split('-')+p2.split('-')
    #out.append('PE'+str(i))
    fo.write('\t'.join(out)+'\n')
    i=i+1

fo.close()
fi.close()

