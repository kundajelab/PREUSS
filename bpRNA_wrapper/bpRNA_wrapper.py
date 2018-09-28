import argparse
import json
import pickle
import pdb
import subprocess
def parse_args():
    parser=argparse.ArgumentParser(description="Wrapper for bpRNA")
    parser.add_argument("--data_json")
    parser.add_argument("--pickle_out")
    parser.add_argument("--text_out",default=None)
    return parser.parse_args()

def get_bpRNA_annotation(input_string):
    #use an intermediate tmp file because bpRNA needs it
    out_tmp=open('tmp.dbn','w')
    out_tmp.write(input_string+'\n')
    out_tmp.close()
    subprocess.call(["bpRNA.pl","tmp.dbn"])
    in_tmp=open("tmp.st",'r')
    annotation=in_tmp.read()
    in_tmp.close()
    return annotation    


def main():
    args=parse_args()
    data=json.load(open(args.data_json,'r'))
    data_dict=dict()
    for item in data['items']:
        #pdb.set_trace()
        cur_id=item['rna_id']
        print("processing structure data for item:"+str(cur_id))
        data_dict[cur_id]=dict()
        sequence_string=item['sequence_string']
        inferred_structure=item['inferred_structure']
        try:
            bootstrap_structures=item['bootstrap_structures']
            data_dict[cur_id]['bootstraps']=dict()
            struct_tally=dict()
            #get a tally of how frequent each bootstrapped structure is
            for struct in bootstrap_structures:
                if struct not in struct_tally:
                    struct_tally[struct]=1
                else:
                    struct_tally[struct]+=1
            bootstrap_count=0
            for struct in struct_tally:
                header=','.join(['>bootstrap',cur_id+'.'+str(bootstrap_count),'frequency:'+str(struct_tally[struct])])
                bpRNA_annotation=get_bpRNA_annotation('\n'.join([header,sequence_string,struct]))
                bootstrap_count+=1
                data_dict[cur_id]['bootstraps'][bpRNA_annotation]=struct_tally[struct]
        except:
            print("no bootstrap structures for "+str(cur_id)+", continuing")
        header=','.join(['>inferred',cur_id])
        bpRNA_annotation=get_bpRNA_annotation('\n'.join([header,sequence_string,inferred_structure]))
        data_dict[cur_id]['inferred']=bpRNA_annotation
    #save to pickle
    with open(args.pickle_out,'wb') as handle:
        pickle.dump(data_dict,handle,protocol=pickle.HIGHEST_PROTOCOL)
    #write to text file
    outf=open(args.text_out,'w')
    for cur_id in data_dict:
        outf.write(">"+cur_id+',inferred\n')
        outf.write(data_dict[cur_id]['inferred']+'\n')
        i=0
        if 'bootstraps' in data_dict[cur_id]: 
            for entry in data_dict[cur_id]['bootstraps']:
                outf.write('>'+cur_id+'.'+str(i)+',bootstrap,count='+str(data_dict[cur_id]['bootstraps'][entry])+'\n')
                outf.write(entry)
                i+=1

        
if __name__=="__main__":
    main()
    
