import _pickle as pickle
import sys
from scipy import stats
import numpy
import sys

########################################
R=0
Connector='.And.'
########################################

print('''$1: Index file\n$2: Output file''')
INDEX_DATA_FILE=sys.argv[1]
OUT_FILE=sys.argv[2]
SPLIT=':'

class Used_Data:
    import numpy
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




print('loading...')
fdata = open(INDEX_DATA_FILE,'rb')
data = pickle.load(fdata)
fdata.close()
print("loading done !")

header=data.header
NETGENE=set()
for edge in data.PAIR_POOL:
    p1=edge.split(SPLIT)[0]
    p2=edge.split(SPLIT)[1]
    NETGENE.add(p1)
    NETGENE.add(p2)
    

def rmOut(X):
    import numpy as np
    X=np.array(X)
    Q3=np.quantile(X,0.75)
    Q1=np.quantile(X,0.25)
    RANGE=Q3-Q1
    UP=Q3+1.5*RANGE
    LW=Q1-1.5*RANGE
    OUT=X[np.argwhere(X<UP)]
    OUT=OUT[np.argwhere(OUT>LW)]
    return(OUT)

def calABCD(X, Y, X_base, Y_base):
    import numpy as np
    X=np.array(X)
    Y=np.array(Y)
    X_base=X_base
    Y_base=Y_base
    X_delta=X-X_base
    Y_delta=Y-Y_base
    A = np.sum(X_delta * Y_delta)
    B = np.sum(X_delta ** 2)
    C = np.sum(Y_delta ** 2)
    ############################
    D = A / np.sqrt( B * C )
    OUT={}
    OUT['A']=A
    OUT['B']=B
    OUT['C']=C
    OUT['D']=D
    return(OUT)



def calILS(X, Y, r=0):
    import numpy as np
    import scipy.stats as st
    X=np.array(X)
    Y=np.array(Y)
    ###########################
    #Ensure positive value
    if np.min(X)<0 or np.min(Y)<0:
        X=X-np.min(rmOut(X))
        Y=Y-np.min(rmOut(Y))
        X[np.argwhere[X<0]]=0
        Y[np.argwhere[Y<0]]=0
    ###########################
    r=r
    N=len(X)
    ############################
    X_mean = np.mean(X)
    Y_mean = np.mean(Y)
    ###########################
    X_base = X_mean * r
    Y_base = Y_mean * r
    X_delta = X - X_base
    Y_delta = Y - Y_base
    ###########################
    ABCD = calABCD(X, Y, X_base, Y_base)
    ###########################
    D = ABCD['D']
    D_plus = ( ABCD['A'] + X_delta * Y_delta ) / np.sqrt( (ABCD['B'] + X_delta**2) * (ABCD['C'] + Y_delta**2) )
    ###########################
    M = D_plus-D
    S = (1-D**2)/(N-1)
    ###########################
    ILS = M / S
    return( ILS )



def calZ(a, b, r):
    ILS = calILS(a, b, r) 
    return(ILS)



def SINGLE(p1, p2, p1_exp, p2_exp):

    this_out=calZ(p1_exp, p2_exp, R)
    Z=[]
    for z in this_out:
        if str(z)=='nan':
            z=0
        Z.append(str(z))
    out=p1+Connector+p2+'\t'+'\t'.join(Z)+'\n'
    return(out)




open(OUT_FILE,'w').write('\t'.join(header)+'\n')

for edge in data.EDGE:
    ps=edge.split(SPLIT)
    p1=ps[0]
    p2=ps[1]
    if p1 in data.EXP and p2 in data.EXP:
        p1_exp=data.EXP[p1]
        p2_exp=data.EXP[p2]
         
        p1_this_list=data.EXP[p1]
        p2_this_list=data.EXP[p2]
        
        this_out=SINGLE(p1, p2, p1_exp,p2_exp)
        open(OUT_FILE,'a').write(this_out)


        











