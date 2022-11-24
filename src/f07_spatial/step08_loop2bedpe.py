


fi=open('ics_out/UNIQ_LOOP.txt')
fo=open('scATAC_RefLoop.bed','w')

fi.readline()
i=1
for line in fi:
    #seq=line.rstrip().split('\t')
    info=line.split('\t')
    p1=info[0]
    p2=info[1]
    out=p1.split('-')+p2.split('-')
    #out.append('PE'+str(i))
    fo.write('\t'.join(out)+'\n')
    i=i+1

fo.close()
fi.close()

