import argparse
import numpy as np
import pandas as pd
import pickle

ltrdict = {'a':[1,0,0,0],
           'c':[0,1,0,0],
           'g':[0,0,1,0],
           'u':[0,0,0,1],
           'n':[0,0,0,0],
           'A':[1,0,0,0],
           'C':[0,1,0,0],
           'G':[0,0,1,0],
           'U':[0,0,0,1],
           'N':[0,0,0,0]}

structdict={'S':[1,0,0,0,0,0],
            'H':[0,1,0,0,0,0],
            'B':[0,0,1,0,0,0],
            'I':[0,0,0,1,0,0],
            'E':[0,0,0,0,1,0],
            'M':[0,0,0,0,0,1],
            'N':[0,0,0,0,0,0]}

    

def parse_args():
    parser=argparse.ArgumentParser(description="aggregate ouptuts from bpRNA")
    parser.add_argument("--neil1_bprna")
    parser.add_argument("--neil1_editing_df")
    parser.add_argument("--ttyh2_bc_bprna")
    parser.add_argument("--ttyh2_bc_editing_df")
    parser.add_argument("--ttyh2_ecs_bprna")
    parser.add_argument("--ttyh2_ecs_editing_df")    
    parser.add_argument("--ajuba_bprna")
    parser.add_argument("--ajuba_editing_df")
    parser.add_argument("--outf")
    return parser.parse_args()

def one_hot_encode_seq(seqs):
    return np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs])

def one_hot_encode_struct(structs):
    return np.array([[structdict.get(x) for x in struct] for struct in structs])
    #return np.array([[structdict.get(x,[0,0,0,0,0]) for x in struct] for struct in structs])

def pad(entries,desired_length):
    padded=[]
    for entry in entries:
        to_add=desired_length-len(entry)
        if to_add%2==0:
            left_flank=to_add//2
            right_flank=to_add//2
        else:
            left_flank=to_add//2
            right_flank=(to_add+1)//2
        new_entry='N'*left_flank+entry+'N'*right_flank
        assert len(new_entry)==desired_length
        padded.append(new_entry)
    return padded 

def process_bpRNA(bpRNA_input,substrate):
    data=[line.split('\n') for line in open(bpRNA_input,'r').read().strip().split('>')]
    max_seq_length=0
    max_struct_length=0 
    names=[]
    seqs=[]
    structs=[] 
    for isoform in data:
        if len(isoform)<7:
            continue
        iso_name=isoform[0]
        if iso_name.endswith('computational')==False:
            continue
        iso_name=substrate+'_'+str(int(iso_name.split(',')[0]))
        names.append(iso_name)
        iso_seq=isoform[4]
        seqs.append(iso_seq)
        if len(iso_seq)> max_seq_length:
            max_seq_length=len(iso_seq) 
        iso_struct=isoform[6]
        structs.append(iso_struct)
        if len(iso_struct)> max_struct_length:
            max_struct_length=len(iso_struct)
    #pad sequences and structs to max
    return names,seqs,structs,max_seq_length,max_struct_length

def make_y_dict(df,substrate):
    cur_dict={}
    for index,row in df.iterrows():
        cur_dict[substrate+"_"+str(int(row['rna_id']))]=row['editing_value']
    return cur_dict 

def main():
    args=parse_args()

    neil1_editing=make_y_dict(pd.read_csv(args.neil1_editing_df)[['rna_id','editing_value']],'NEIL1')
    ttyh2_bc_editing=make_y_dict(pd.read_csv(args.ttyh2_bc_editing_df)[['rna_id','editing_value']],'TTYH2_BC')
    ttyh2_ecs_editing=make_y_dict(pd.read_csv(args.ttyh2_ecs_editing_df)[['rna_id','editing_value']],'TTYH2_ECS')
    ajuba_editing=make_y_dict(pd.read_csv(args.ajuba_editing_df)[['rna_id','editing_value']],'AJUBA')
    y_labels=neil1_editing
    y_labels.update(ttyh2_bc_editing)
    y_labels.update(ttyh2_ecs_editing)
    y_labels.update(ajuba_editing)
    print("made y labels")

    neil1_names,neil1_seq,neil1_struct,max_len_neil1_seq,max_len_neil1_struct=process_bpRNA(args.neil1_bprna,'NEIL1')
    ttyh2_bc_names,ttyh2_bc_seq,ttyh2_bc_struct,max_len_ttyh2_bc_seq,max_len_ttyh2_bc_struct=process_bpRNA(args.ttyh2_bc_bprna,'TTYH2_BC')
    ttyh2_ecs_names,ttyh2_ecs_seq,ttyh2_ecs_struct,max_len_ttyh2_ecs_seq,max_len_ttyh2_struct=process_bpRNA(args.ttyh2_ecs_bprna,'TTYH2_ECS')
    ajuba_names,ajuba_seq,ajuba_struct,max_len_ajuba_seq,max_len_ajuba_struct=process_bpRNA(args.ajuba_bprna,'AJUBA')

    max_seq_length=max([max_len_neil1_seq, max_len_ttyh2_bc_seq, max_len_ttyh2_ecs_seq, max_len_ajuba_seq ])
    max_struct_length=max([max_len_neil1_struct, max_len_ttyh2_bc_struct, max_len_ttyh2_ecs_seq, max_len_ajuba_struct])
    
    print("max seq length:"+str(max_seq_length))
    print("max struct length:"+str(max_struct_length))
    neil1_seq=one_hot_encode_seq(pad(neil1_seq,max_seq_length))
    ttyh2_bc_seq=one_hot_encode_seq(pad(ttyh2_bc_seq,max_seq_length))
    ttyh2_ecs_seq=one_hot_encode_seq(pad(ttyh2_ecs_seq,max_seq_length))
    ajuba_seq=one_hot_encode_seq(pad(ajuba_seq,max_seq_length))
    print("one-hot-encoded sequences") 
    neil1_struct=one_hot_encode_struct(pad(neil1_struct,max_struct_length))
    print(str(neil1_struct.shape))
    ttyh2_bc_struct=one_hot_encode_struct(pad(ttyh2_bc_struct,max_struct_length))
    print(str(ttyh2_bc_struct.shape))
    ttyh2_ecs_struct=one_hot_encode_struct(pad(ttyh2_ecs_struct,max_struct_length))
    print(str(ttyh2_ecs_struct.shape))
    ajuba_struct=one_hot_encode_struct(pad(ajuba_struct,max_struct_length))
    print(str(ajuba_struct.shape))
    print("one-hot-encoded structs") 
    X_vals={}
    for i in range(len(neil1_names)):
        X_vals[neil1_names[i]]=[neil1_seq[i],neil1_struct[i]]
    for i in range(len(ttyh2_bc_names)):
        X_vals[ttyh2_bc_names[i]]=[ttyh2_bc_seq[i],ttyh2_bc_struct[i]]
    for i in range(len(ttyh2_ecs_names)):
        X_vals[ttyh2_ecs_names[i]]=[ttyh2_ecs_seq[i],ttyh2_ecs_struct[i]]
    for i in range(len(ajuba_names)):
        X_vals[ajuba_names[i]]=[ajuba_seq[i],ajuba_struct[i]]

    #save X & y vals to a pickle
    entries={}
    entries['X']=X_vals
    entries['y']=y_labels
    pickle.dump(entries,open(args.outf,'wb'))
    
    
                
    
if __name__=="__main__":
    main()
    
