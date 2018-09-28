import pandas as pd
import pickle
import pdb
import json 
import argparse

def parse_args():
    parser=argparse.ArgumentParser(description="generate feature matrix for adar edited RNA")
    parser.add_argument("--rna_lib_structure_summary_json")
    parser.add_argument("--bpRNA_pickle")
    parser.add_argument("--outf")
    parser.add_argument("--annotate_bootstraps",action='store_true',default=False)
    return parser.parse_args()

def format_id(rna_id):
    if rna_id<10:
        rna_id="00"+str(rna_id)
    elif rna_id<100:
        rna_id="0"+str(rna_id)
    else:
        rna_id=str(rna_id) 
    return rna_id

def get_bprna_feature_labels(pos,features,editing_site):
    pos=int(pos)
    editing_site=int(editing_site)
    try:
        cur_feat=features[pos]
        prev_feat=None
        next_feat=None
        #previous feature
        for i in range(pos-1,0,-1):
            if features[i]!=cur_feat:
                prev_feat=features[i]
                break
        #subsequent feature
        for i in range(pos+1,len(features)):
            if features[i]!=cur_feat:
                next_feat=features[i]
                break
        pos_to_edit=features[min([pos,editing_site]):max([pos,editing_site])+1]
        if (len(set(pos_to_edit))==1):
            cur_feat_same_as_edit=1
        else:
            cur_feat_same_as_edit=0
    except:
        cur_feat=None
        prev_feat=None
        next_feat=None
        cur_feat_same_as_edit=None
    return cur_feat,prev_feat,next_feat,cur_feat_same_as_edit

def get_structure_length(features,editing_site,structure_type):
    structure_start=None
    structure_end=None
    for i in range(editing_site,0,-1):
        if features[i]==structure_type:
            structure_end=i
            break
    for i in range(structure_end,0,-1):
        if features[i]!=structure_type:
            structure_start=i
            break
    structure_length=structure_end-structure_start
    return structure_length

def get_downstream_nonstem_info(features,editing_site,inferred):
    #find the first non-stem feature 3' of editing site (might include the editing site)
    feat_3prime_of_edit_site=None
    hook_pos=None
    for pos in range(editing_site,len(features)):
        if features[pos]!="S":
            feat_3prime_of_edit_site=features[pos]
            hook_pos=pos
            break
    #get the start & end position of this feature
    feat_start=None
    feat_end=None
    for i in range(hook_pos,0,-1):
        if features[i]!=feat_3prime_of_edit_site:
            feat_start=i+1
            break
    for i in range(hook_pos,len(features)):
        if features[i]!=feat_3prime_of_edit_site:
            feat_end=i-1
            break
    try:
        feat_3prime_of_edit_site_length=feat_end-feat_start+1
    except:
        return None,None,None,None,None,None
    feat_3prime_of_edit_site_length_5prime=0

        
    #now, we know the type of feature, the start/end positions
    #iterate through bpRNA annotation to get additional info (i.e. closing pairs)
    for entry in inferred[7::]:
        if entry.startswith(feat_3prime_of_edit_site):
            entry_tokens=entry.strip().split(' ')
            entry_pos=[int(i) for i in entry_tokens[1].split('..')]
            #Careful! comparing 0-indexed positions to 1-indexing 
            if (((feat_start+1)==entry_pos[0]) and ((feat_end+1)==entry_pos[1])):
                #this is the feature we want
                closing_pair1=entry_tokens[4]
                if len(entry_tokens)>5:
                    closing_pair2=entry_tokens[-1]
                    to_find=None
                else:
                    #find the complementary entry for the internal loop
                    try:
                        tofind_base=entry_tokens[0].split('.')[0]
                        found_suffix=entry_tokens[0].split('.')[1]
                        if found_suffix=="1":
                            to_find=tofind_base+'.2'
                        else:
                            to_find=tofind_base+'.1'
                    except:
                        to_find=None
                        closing_pair2=None
                break
    if(to_find!=None):
        for entry in inferred[7::]:
            if entry.startswith(to_find):
                closing_pair2=entry.strip().split(' ')[-1]
                pos_info=entry.strip().split(' ')[1]
                pos_start=pos_info.split('..')[0]
                pos_end=pos_info.split('..')[1]
                feat_3prime_of_edit_site_length_5prime=int(pos_end)-int(pos_start)+1 
                break
    return feat_3prime_of_edit_site,feat_3prime_of_edit_site_length,feat_3prime_of_edit_site_length_5prime,closing_pair1,closing_pair2,feat_end

def annotate_structure(editing_levels,bprna_data):
    structure_features=dict() 
    for cur_id in bprna_data:
        inferred=bprna_data[cur_id]['inferred'].split('\n')
        features=inferred[5]
        cur_sequence=inferred[3]
        cur_structure=inferred[4]
        
        #Annotate Editing Site Features
        #IMPORTANT -- DATA IS PROVIDED 1-INDEXED. SHIFT TO 0-INDEX FOR USE IN INDEXING FEATURE ARRAY
        try:
            editing_site=int(editing_levels[cur_id]['site'])-1
        except:
            pdb.set_trace() 
        editing_feature=features[editing_site]
        mfeat1=None
        mfeat1_prev=None
        mfeat1_next=None
        mfeat1_same_as_edit=None
        mfeat2=None
        mfeat2_prev=None
        mfeat2_next=None
        mfeat2_same_as_edit=None
        
        #Annotate first & second  mutation
        #SHIFT mp1 & mp2 to 0-indexing as well
        mp1=editing_levels[cur_id]['mut_info']['mp1']
        mp2=editing_levels[cur_id]['mut_info']['mp2']
        if (mp1!=None):
            mp1=int(mp1)-1
            mfeat1,mfeat1_prev,mfeat1_next,mfeat1_same_as_edit=get_bprna_feature_labels(mp1,features,editing_site)
        if (mp2!=None):
            mp2=int(mp2)-1
            mfeat2,mfeat2_prev,mfeat2_next,mfeat2_same_as_edit=get_bprna_feature_labels(mp2,features,editing_site)
            
        #get the stem length
        stem_length=get_structure_length(features,editing_site,'S') 
        #get the hairpin length
        hairpin_length=get_structure_length(features,editing_site,'H')
        #length of internal loop downstream of editing site, type of loop, closing pair 
        feat_3prime_of_edit_site,feat_3prime_of_edit_site_length,feat_3prime_of_edit_site_length_5prime,closing_pair1,closing_pair2,feat_3prime_of_edit_site_endpos=get_downstream_nonstem_info(features,editing_site,inferred)
        feat_3prime_of_edit_site_distal,feat_3prime_of_edit_site_length_distal,feat_3prime_of_edit_site_length_5prime_distal,closing_pair1_distal,closing_pair2_distal,endpos_distal=get_downstream_nonstem_info(features,feat_3prime_of_edit_site_endpos,inferred)
        #store all to dict
        structure_features[cur_id]=dict()
        structure_features[cur_id]['editing_feature']=editing_feature
        structure_features[cur_id]['mfeat1']=mfeat1
        structure_features[cur_id]['mfeat1_prev']=mfeat1_prev
        structure_features[cur_id]['mfeat1_next']=mfeat1_next
        structure_features[cur_id]['mfeat1_same_as_edit']=mfeat1_same_as_edit
        structure_features[cur_id]['mfeat2']=mfeat2
        structure_features[cur_id]['mfeat2_prev']=mfeat2_prev
        structure_features[cur_id]['mfeat2_next']=mfeat2_next
        structure_features[cur_id]['mfeat2_same_as_edit']=mfeat2_same_as_edit
        structure_features[cur_id]['stem_length']=stem_length
        structure_features[cur_id]['hairpin_length']=hairpin_length
        structure_features[cur_id]['feat_3prime_e']=feat_3prime_of_edit_site
        structure_features[cur_id]['feat_3prime_e_length']=feat_3prime_of_edit_site_length
        structure_features[cur_id]['feat_3prime_e_length_5prime']=feat_3prime_of_edit_site_length_5prime
        structure_features[cur_id]['feat_3prime_e_cp1']=closing_pair1
        structure_features[cur_id]['feat_3prime_e_cp2']=closing_pair2
        structure_features[cur_id]['feat_3prime_e_distal']=feat_3prime_of_edit_site_distal
        structure_features[cur_id]['feat_3prime_e_length_distal']=feat_3prime_of_edit_site_length_distal
        structure_features[cur_id]['feat_3prime_e_length_5prime_distal']=feat_3prime_of_edit_site_length_5prime_distal
        structure_features[cur_id]['feat_3prime_e_cp1_distal']=closing_pair1_distal
        structure_features[cur_id]['feat_3prime_e_cp2_distal']=closing_pair2_distal

    return structure_features 
    
def write_feature_matrix(editing_levels,structure_dict,outf):
    outf=open(outf,'w')
    header=None
    for cur_id in editing_levels.keys():
        mut_info_keys=list(editing_levels[cur_id]['mut_info'].keys())
        struct_info_keys=list(structure_dict[cur_id].keys())
        editing_level=editing_levels[cur_id]['level']
        if header==None:
            header=['cur_id','editing_level']+mut_info_keys+struct_info_keys
            header='\t'.join([str(i) for i in header])
            outf.write(header+'\n')
        mut_info=editing_levels[cur_id]['mut_info']
        struct_info=structure_dict[cur_id]
        outf.write(str(cur_id)+'\t'+str(editing_level))
        for keyname in mut_info_keys:
            outf.write('\t'+str(mut_info[keyname]))
        for keyname in struct_info_keys:
            outf.write('\t'+str(struct_info[keyname]))
        outf.write('\n')
                           
def get_editing_info(editing_levels_file):
    data=json.load(open(editing_levels_file,'r'))
    editing_levels_dict=dict()
    for item in data['items']:
        cur_id=item['rna_id']
        editing_site=int(item['A-to-I_editing_site'])
        ave_editing_level=item['A-to-I_editing_level']
        mutation_syntax=item['mutation_syntax']
        if mutation_syntax==None:
            mutation_syntax="WT"
        if mutation_syntax=="None":
            mutation_syntax="WT"
        editing_levels_dict[cur_id]=dict()
        editing_levels_dict[cur_id]['level']=ave_editing_level
        editing_levels_dict[cur_id]['site']=editing_site
        editing_levels_dict[cur_id]['mutation_syntax']=mutation_syntax 
    return editing_levels_dict

def annotate_mutation(mut,editing_site):
    mp1=mut[:-4]
    adist1=editing_site-int(mp1)
    mref1=mut[-4]
    malt1=mut[-1]
    mtype1="mismatch"
    return mp1,adist1,mref1,malt1,mtype1

def get_mut_info(editing_levels_dict):
    for cur_id in editing_levels_dict:
        mut_syntax=editing_levels_dict[cur_id]['mutation_syntax'].split(',')
        editing_site=editing_levels_dict[cur_id]['site']
        #keep track of mutation position information
        mp1=None
        mp2=None
        mref1=None
        mref2=None
        malt1=None
        malt2=None
        mtype1=None
        mtype2=None
        adist1=None
        adist2=None

        if len(mut_syntax)==2:
            #annotate the first mutation
            mp1,adist1,mref1,malt1,mtype1=annotate_mutation(mut_syntax[0],editing_site)
            mp2,adist2,mref2,malt2,mtype2=annotate_mutation(mut_syntax[1],editing_site)
        elif "to" in mut_syntax[0]:
            #This is a single point mutation
            mp1,adist1,mref1,malt1,mtype1=annotate_mutation(mut_syntax[0],editing_site)            
        elif "Indel" in mut_syntax[0]:
            #This is an indel 
            mtype1="indel"
            mp1=mut_syntax[0][5::].split('-')[0]
        elif "WT" in mut_syntax[0]:
            #This is wild type
            mtype1="wt"
        else:
            pdb.set_trace()
        editing_levels_dict[cur_id]['mut_info']={"mp1":mp1,
                               "mp2":mp2,
                               "mref1":mref1,
                               "mref2":mref2,
                               "malt1":malt1,
                               "malt2":malt2,
                               "mtype1":mtype1,
                               "mtype2":mtype2,
                               "adist1":adist1,
                               "adist2":adist2}
    return editing_levels_dict

def inferred_base_in_stem_freq(editing_levels_dict,bprna_data,outf,feat_type):
    entries=editing_levels_dict.keys()
    stem_freq=dict()
    for entry in entries:
        bootstraps=bprna_data[entry]['bootstraps']        
        struct_len=len(bprna_data[entry]['inferred'].split('\n')[5])
        stem_count=[0]*struct_len
        num_bootstraps=0
        #weight frequencies by the number of times the bootstrap is observed. 
        for bootstrap in bootstraps:
            bootstrap_count=bootstraps[bootstrap]
            num_bootstraps+=bootstrap_count
            structure=bootstrap.split('\n')[5]
            for base in range(struct_len):
                if structure[base]==feat_type:
                    stem_count[base]+=bootstrap_count
        stem_freq[entry]=[float(i)/num_bootstraps for i in stem_count]
    #write the output file
    outf=open(outf+'.'+feat_type+'.freq.txt','w')
    outf.write('RNA_ID\tEditingLevel\t'+'\t'.join([feat_type+'.base'+str(i) for i in range(struct_len)])+'\n')
    for entry in editing_levels_dict:
        outf.write(entry+'\t'+str(editing_levels_dict[entry]['level'])+'\t'+'\t'.join([str(i) for i in stem_freq[entry]])+'\n')


def main():
    args=parse_args()
    editing_levels_file=args.rna_lib_structure_summary_json
    bprna_data=pickle.load(open(args.bpRNA_pickle,'rb'))
    outf=args.outf
    
    editing_levels_dict=get_editing_info(editing_levels_file)
    editing_levels_dict=get_mut_info(editing_levels_dict)
    
    #what fraction of bases in bootstrapped samples are in specified_feature?
    if (args.annotate_bootstraps==True):
        inferred_base_in_stem_freq(editing_levels_dict,bprna_data,outf,'S')
        inferred_base_in_stem_freq(editing_levels_dict,bprna_data,outf,'I')
        inferred_base_in_stem_freq(editing_levels_dict,bprna_data,outf,'B')
        inferred_base_in_stem_freq(editing_levels_dict,bprna_data,outf,'H')
    

    #annotate inferred structure
    structure_dict=annotate_structure(editing_levels_dict,bprna_data)
    write_feature_matrix(editing_levels_dict,structure_dict,outf)
    
if __name__=="__main__":
    main()
