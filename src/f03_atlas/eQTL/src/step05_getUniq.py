
import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')


fa=open('/home/disk/database/data/ICS_scHIC/f03_atlas/eQTL/LIST.txt')
i=1
old=set()
for line in fa:
    seq=line.rstrip().split('\t')
    input_path='/home/disk/database/data/ICS_scHIC/f03_atlas/eQTL/GTEx_Analysis_v8_eQTL/'+seq[0]+'.format/tag.bed'
    if input_path != sys.argv[1]:
        print(input_path)
        this_fi=open(input_path)    
        for lll in this_fi:
            old.add(lll.rstrip())

#print(len(old))

for line in fi:
    if line.rstrip() not in old:
        fo.write(line)

fo.close()



