import sys
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.longrange','w')
for line in fi:
    seq=line.rstrip().split('\t')
    this_out=seq[0]+'\t'+seq[1]+'\t'+seq[2]+'\t'+seq[3]+':'+seq[4]+'-'+seq[5]+','+seq[6]+'\n'
    fo.write(this_out)
fi.close()
fo.close()


