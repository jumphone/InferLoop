import sys

fi=open(sys.argv[1])
fo1=open(sys.argv[2],'w')
fo2=open(sys.argv[3],'w')
SLOP=int(sys.argv[4])

BED=set()
for line in fi:
    seq=line.rstrip().split('\t') 
    new_start_1=seq[1] #str(max([int(seq[1])-SLOP,0]))
    new_end_1=seq[2] #str(int(seq[2])+SLOP)
    new_start_2=str(max([int(seq[4])-SLOP,0]))
    new_end_2=str(int(seq[5])+SLOP)
    
    this_snp=seq[0]+'x'+seq[1]+'x'+seq[2]
    this_gene=seq[6]

    this_bed1=seq[0]+'\t'+new_start_1+'\t'+new_end_1+'\t'+this_snp
    this_bed2=seq[3]+'\t'+new_start_2+'\t'+new_end_2+'\t'+this_gene
    BED.add(this_bed1)
    BED.add(this_bed2)
    this_loop1=seq[0]+'-'+new_start_1+'-'+new_end_1
    this_loop2=seq[3]+'-'+new_start_2+'-'+new_end_2
    fo2.write(this_loop1+'\t'+this_loop2+'\t'+this_snp+'\t'+this_gene+'\t'+str(abs(int(seq[5])-int(seq[2])))+'\n')

OUT=[]
for one in BED:
   one=one.split('\t')
   OUT.append([one[0],int(one[1]),int(one[2]),one[3]])

OUT.sort()

for one in OUT:
    fo1.write(one[0]+'\t'+str(one[1])+'\t'+str(one[2])+'\t'+one[3]+'\n')

fo1.close()
fo2.close()
