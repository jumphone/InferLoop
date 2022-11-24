
import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')
OLD=set()
for line in fi:
    seq=line.rstrip().split('\t')
    this_net=[seq[0],seq[1]]
    this_net.sort()
    this_tmp = ':'.join(this_net)
    if this_tmp not in OLD:
        fo.write(line)
        OLD.add(this_tmp)

fo.close()
fi.close()
