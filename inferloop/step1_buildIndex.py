import sys
import _pickle as pickle
import numpy

print('\n$1 Net_file (col1: ID; col2: ID)\n$2 Signal_matrix (col: bins; row: ID)\n$3 Output_index_file\n')

NETWORK_PATH=sys.argv[1]
EXP_POOL_PATH=sys.argv[2]
OUT_DATA_PATH=sys.argv[3]
SPLIT=':'

class Used_Data:
    def __init__(self, EXP, EXP_LENGTH, EDGE, header):
        self.EXP = EXP
        self.EXP_LENGTH = EXP_LENGTH
        self.header = header
        ##############
        self.EDGE = EDGE


EDGE=[]
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
        EDGE.append(edge)
        POINT.add(p1)
        POINT.add(p2)
fnet.close()

EXP={}
EXP_GENE=set()
EXP_LENGTH=0
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
            EXP_LENGTH=len(tmp)
            EXP[seq[0]]=tmp
            if len(header)==len(EXP[seq[0]])+1:
                header=header[1:]

fpool.close()
print('ncol(mat)='+str(EXP_LENGTH))

NETGENE=set()
for edge in EDGE:
    ps=edge.split(SPLIT)
    p1=ps[0]
    p2=ps[1]
    if 1==1:
        p1_exp=EXP[p1]
        p2_exp=EXP[p2]
        NETGENE.add(p1)
        NETGENE.add(p2)

OUTEXP={}
for gene in EXP:
    if gene in NETGENE:
        OUTEXP[gene]=EXP[gene]

print('nrow(net)='+str(len(EDGE)))
print('\n')
fo=open(OUT_DATA_PATH,'wb')
data=Used_Data( OUTEXP, EXP_LENGTH, EDGE, header)
pickle.dump(data,fo)
fo.close()



