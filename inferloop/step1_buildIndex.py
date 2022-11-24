import sys
import _pickle as pickle
import numpy

print('\n$1 Net_file (col1: ID; col2: ID)\n$2 Signal_matrix (col: bins; row: ID)\n$3 Output_index_file\n')

NETWORK_PATH=sys.argv[1]
EXP_POOL_PATH=sys.argv[2]
OUT_DATA_PATH=sys.argv[3]
SPLIT=':'

class Used_Data:
    def __init__(self, EXP, POOL_LENGTH, PAIR_POOL, header):
        self.EXP = EXP
        self.POOL_LENGTH = POOL_LENGTH
        self.PAIR_POOL = PAIR_POOL
        self.header = header
        ##############
        self.EDGE = []
        for this_edge in self.PAIR_POOL:
            self.EDGE.append(this_edge)
        self.EDGE.sort()
        ##############
        self.EXP_MEAN = {}
        for gene in EXP:
            self.EXP_MEAN[gene]=numpy.mean(EXP[gene])
        ##############
        self.EXP_STD = {}
        for gene in EXP:
            self.EXP_STD[gene]=numpy.std(EXP[gene])
        ##############
        self.EXP_MIN = {}
        for gene in EXP:
            self.EXP_MIN[gene]=numpy.min(EXP[gene])
        ##############
        self.EXP_MAX = {}
        for gene in EXP:
            self.EXP_MAX[gene]=numpy.max(EXP[gene])


EDGE=set()
POINT=set()
fnet=open(NETWORK_PATH)
for line in fnet:
    seq=line.rstrip().split('\t')
    if line[0]!='#':
        p1=seq[0]
        p2=seq[1]
        ######################
        edge=[p1,p2]
        edge=SPLIT.join(edge)
        ########################
        EDGE.add(edge)
        POINT.add(p1)
        POINT.add(p2)
fnet.close()

EXP={}
EXP_GENE=set()
POOL_LENGTH=0
fpool=open(EXP_POOL_PATH)
header=fpool.readline().rstrip().split('\t')
for line in fpool:
        seq=line.rstrip().split('\t')
        if seq[0] in POINT:
            EXP_GENE.add(seq[0])
            tmp=[]
            for one in seq[1:]:
                try:
                    tmp.append(float(one))
                except IOError:
                    print(one)
                    tmp.append(0.0)
            POOL_LENGTH=len(tmp)
            EXP[seq[0]]=tmp
            if len(header)==len(EXP[seq[0]])+1:
                header=header[1:]

fpool.close()
print(POOL_LENGTH)



EDGE=set()
POINT=set()
fnet=open(NETWORK_PATH)
for line in fnet:
    seq=line.rstrip().split('\t')
    if line[0]!='#':
        p1=seq[0]
        p2=seq[1]
        if p1 in EXP_GENE and p2 in EXP_GENE:
            ######################
            edge=[p1,p2]
            edge=SPLIT.join(edge)
            ########################
            EDGE.add(edge)
            POINT.add(p1)
            POINT.add(p2)

fnet.close()



from scipy import stats
PAIR_POOL={}
NETGENE=set()
for edge in EDGE:
    ps=edge.split(SPLIT)
    p1=ps[0]
    p2=ps[1]
    if 1==1:
        p1_exp=EXP[p1]
        p2_exp=EXP[p2]
        PAIR_POOL[edge]=[float(len(p1_exp))]
        NETGENE.add(p1)
        NETGENE.add(p2)

OUTEXP={}
for gene in EXP:
    if gene in NETGENE:
        OUTEXP[gene]=EXP[gene]

print(len(PAIR_POOL))
fo=open(OUT_DATA_PATH,'wb')
data=Used_Data( OUTEXP, POOL_LENGTH, PAIR_POOL, header)
pickle.dump(data,fo)
fo.close()



