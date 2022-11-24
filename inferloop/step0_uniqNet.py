


print('\n$1 Net_file (col1: ID; col2: ID)\n$2 Output_net_file \n')


import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')

SPLIT=':'

OLD=set()
for line in fi:
    seq=line.rstrip().split('\t')
    this_net=[seq[0],seq[1]]
    this_net.sort()
    this_tmp = SPLIT.join(this_net)
    if this_tmp not in OLD:
        fo.write(line)
        OLD.add(this_tmp)

fo.close()
fi.close()
