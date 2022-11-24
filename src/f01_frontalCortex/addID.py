import sys

fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')

i=1
for line in fi:
    seq=line.rstrip().split('\t')
    this_id='ID'+str(i)+'E'
    this_out=seq[0:3]+[this_id]+seq[3:]
    out='\t'.join(this_out)+'\n'
    fo.write(out)
    i=i+1
fo.close()
fi.close()
