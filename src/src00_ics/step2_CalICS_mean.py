import _pickle as pickle
import sys
from scipy import stats
import numpy
import multiprocessing

print('''$1: Index Data\n$2: Output file''')

INDEX_DATA_FILE=sys.argv[1]
OUT_FILE=sys.argv[2]
SPLIT=':'

class Person_PCC_Data:
    import numpy
    def __init__(self, EXP, POOL_LENGTH, PCC_POOL,header):
        self.EXP = EXP
        self.POOL_LENGTH = POOL_LENGTH
        self.PCC_POOL = PCC_POOL
        self.header = header
        ##############
        self.EDGE = []
        for this_edge in self.PCC_POOL:
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
for edge in data.PCC_POOL:
    p1=edge.split(SPLIT)[0]
    p2=edge.split(SPLIT)[1]
    NETGENE.add(p1)
    NETGENE.add(p2)
    



from sklearn.preprocessing import StandardScaler


def SINGLE(p1, p2,p1_old_exp,p2_old_exp, p1_this_list, p2_this_list, p1_dict, p2_dict, pcc_old,pcc_length):
    if 1==1:        
        p1_this_list = StandardScaler(with_mean=False).fit_transform(numpy.array(p1_this_list).reshape(-1,1))[:,0]
        p2_this_list = StandardScaler(with_mean=False).fit_transform(numpy.array(p2_this_list).reshape(-1,1))[:,0]
        
        Z=[]
        j=0
        while j<len(header):
            p1_this=p1_this_list[j]
            p2_this=p2_this_list[j]
            if 1==1:
                z = (p1_this + p2_this) / 2
            Z.append(z)
            j+=1

        Z = StandardScaler().fit_transform(numpy.array(Z).reshape(-1,1))[:,0]


        FZ=[]
        for z in Z:
            if str(z)=='nan':
                z=0
            FZ.append(str(z))
        out=p1+'.And.'+p2+'\t'+'\t'.join(FZ)+'\n'
        return(out)




open(OUT_FILE,'w').write('\t'.join(header)+'\n')

for edge in data.EDGE:
    ps=edge.split(SPLIT)
    p1=ps[0]
    p2=ps[1]
    if p1 in data.EXP and p2 in data.EXP:
        p1_old_exp=data.EXP[p1]
        p2_old_exp=data.EXP[p2]
         
        p1_this_list=data.EXP[p1]
        p2_this_list=data.EXP[p2]
        
        #######################################
        p1_dict={}
        p2_dict={}
        p1_dict['mean']=data.EXP_MEAN[p1]
        p2_dict['mean']=data.EXP_MEAN[p2]        
        p1_dict['std']=data.EXP_STD[p1]
        p2_dict['std']=data.EXP_STD[p2]
        p1_dict['min']=data.EXP_MIN[p1]
        p2_dict['min']=data.EXP_MIN[p2]
        p1_dict['max']=data.EXP_MAX[p1]
        p2_dict['max']=data.EXP_MAX[p2]
        #########################################

        pcc_old = data.PCC_POOL[edge][0]
        pcc_length = data.PCC_POOL[edge][1]
        this_out=SINGLE(p1,p2,p1_old_exp,p2_old_exp,p1_this_list, p2_this_list, p1_dict, p2_dict, pcc_old, pcc_length)
        open(OUT_FILE,'a').write(this_out)


        











