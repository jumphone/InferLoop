import sys
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.ucscInter.bed','w')

i=1
for line in fi:
    seq=line.rstrip().split('\t')
    c1=seq[0]
    s1=seq[1]
    e1=seq[2]
    c2=seq[3]
    s2=seq[4]
    e2=seq[5]
   
    if c1==c2:
        signal=str(min([int(seq[6]),1000]))
        this_chr=c1
        this_start=str(min([int(s1),int(s2),int(e1),int(e2)]))
        this_end=str(max([int(s1),int(s2),int(e1),int(e2)]))
        this_name='inter'+str(i)
        this_s_name='s'+str(i)
        this_t_name='t'+str(i)
        this_col='#7A67EE'
        this_out=[this_chr,this_start,this_end,this_name,signal,'0','.',this_col,c1,s1,e1,
                 this_s_name, '+', c2, s2,e2,this_t_name, '+']
        fo.write('\t'.join(this_out)+'\n')
        i=i+1 


