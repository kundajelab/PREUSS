import sys
#change comma to tab
#convert non-unique feature names to unique feature names 
data=open(sys.argv[1],'r').read().strip().split('\n')
outf=open(sys.argv[1]+'.cleaned','w')
rna_id_dict=dict()
outf.write(data[0].replace(',','\t')+'\n')
for line in data[1::]:
    tokens=line.split(',')
    cur_id=tokens[0]
    if cur_id not in rna_id_dict:
        rna_id_dict[cur_id]=1
    else:
        rna_id_dict[cur_id]+=1
    cur_id='_'.join([cur_id,str(rna_id_dict[cur_id])])
    tokens[0]=cur_id
    outf.write('\t'.join(tokens)+'\n')
    
