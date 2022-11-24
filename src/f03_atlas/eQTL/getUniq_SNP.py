import sys

fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')
V={}
for line in fi:
    seq=line.rstrip().split('\t')
    #print(seq)
    this_tag=seq[0]+'\t'+seq[1]+'\t'+seq[2]
    if this_tag not in V:
        V[this_tag]=[[float(seq[-1]),line]]
    else:
        V[this_tag].append([float(seq[-1]),line])
fi.close()

for one in V:
    tmp=V[one]
    tmp.sort()
    fo.write(tmp[0][1])    
    
fo.close()


