
import sys

fi=open('/home/database/annotation/hg38/gencode.v39.annotation.gtf.pc.sorted.bed')
E2G={}
for line in fi:
    seq=line.rstrip().split('\t') 
    E2G[seq[3].split('.')[0]]=seq[4]
fi.close()
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')


this_tissue=sys.argv[1].split('/')[-1].split('.v8.')[0]

fi.readline()
for line in fi:
    #chr1_845402_A_G_b38     ENSG00000225972.1
    seq=line.rstrip().split('\t')
    this_dist=seq[2]
    #CUT=100000
    if seq[1].split('.')[0] in E2G: #and abs(int(this_dist))<CUT:
        this_gene=E2G[seq[1].split('.')[0]]
        this_info=seq[0].split('_')
        this_pv=seq[6]
        this_slope=seq[7]
        this_out=[this_info[0],str(int(this_info[1])-1),this_info[1],this_gene+'\t'+this_tissue,this_dist,this_pv,this_slope]
        fo.write('\t'.join(this_out)+'\n')


fo.close()


