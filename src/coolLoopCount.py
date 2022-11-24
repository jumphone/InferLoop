import cooler,sys

COOL_FILE=sys.argv[1]
BED_FILE=sys.argv[2]
OUT_FILE=sys.argv[3]

c = cooler.Cooler(COOL_FILE)

fi=open(BED_FILE)
fo=open(OUT_FILE,'w')
for line in fi:
    if line[:3]=='chr':
        seq=line.rstrip().split('\t')
        this_c1=seq[0]
        this_s1=seq[1]
        this_e1=seq[2]
        this_c2=seq[3]
        this_s2=seq[4]
        this_e2=seq[5]
        this_count=sum(sum(c.matrix(balance=False).fetch(this_c1+":"+this_s1+"-"+this_e1,this_c2+":"+this_s2+"-"+this_e2)))
        fo.write('\t'.join(seq[:6])+'\t'+str(this_count)+'\n')
fo.close()
fi.close()





