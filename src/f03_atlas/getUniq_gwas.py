import sys

fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')
V={}
for line in fi:
    seq=line.rstrip().split('\t')
    #print(seq)
    this_tag='\t'.join(seq[:9])
    if this_tag not in V:
        V[this_tag]=[[float(seq[-1]),seq[-2]]]
    else:
        V[this_tag].append([float(seq[-1]),seq[-2]])
fi.close()

for one in V:
    tmp=V[one]
    tmp.sort()
    fo.write(one+'\t'+tmp[0][1]+'\t'+str(tmp[0][0])+'\n')    
    
fo.close()


