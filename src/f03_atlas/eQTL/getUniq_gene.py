fi=open('gencode.v39.annotation.gtf.pc.sorted.changed.bed')
fo=open('gencode.v39.annotation.gtf.pc.sorted.changed.bed.uniqGene.bed','w')
GENE={}
for line in fi:
    seq=line.rstrip().split('\t')
    this_gene=seq[3]
    this_len=int(seq[4])
    if this_gene not in GENE:
       GENE[this_gene]=[[this_len,line]]
    else:
       GENE[this_gene].append([this_len,line])


for this_gene in GENE:
    this_out=GENE[this_gene]
    this_out.sort(reverse=True)
    fo.write(this_out[0][1])

fo.close()
