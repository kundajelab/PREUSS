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
    parser.add_argument("--approach",default="computational",choices=["computational","experimental"])
    parser.add_argument("--source",default="NA") 
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

def get_downstream_nonstem_info(features,editing_site,annotation):
    #find the first non-stem feature 3' of editing site (might include the editing site)
    x1feat_downstream_of_edit_site=None
    hook_pos=None
    try:
        for pos in range(editing_site,len(features)):
            if features[pos]!="S":
                x1feat_downstream_of_edit_site=features[pos]
                hook_pos=pos
                break
    except:
        return None,None,None,None,None,None
    #get the start & end position of this feature
    feat_start=None
    feat_end=None
    for i in range(hook_pos,0,-1):
        if features[i]!=x1feat_downstream_of_edit_site:
            feat_start=i+1
            break
    for i in range(hook_pos,len(features)):
        if features[i]!=x1feat_downstream_of_edit_site:
            feat_end=i-1
            break
    try:
        x1feat_downstream_of_edit_site_length_editing_strand=feat_end-feat_start+1
    except:
        return None,None,None,None,None,None
    x1feat_downstream_of_edit_site_length_complementary_strand=0

        
    #now, we know the type of feature, the start/end positions
    #iterate through bpRNA annotation to get additional info (i.e. closing pairs)
    for entry in annotation[7::]:
        if entry.startswith(x1feat_downstream_of_edit_site):
            entry_tokens=entry.strip().split(' ')
            entry_pos=[int(i) for i in entry_tokens[1].split('..')]
            #Careful! comparing 0-indexed positions to 1-indexing 
            if (((feat_start+1)==entry_pos[0]) and ((feat_end+1)==entry_pos[1])):
                #this is the feature we want
                x1feat_downstream_fo_edit_site_5prime_cp=entry_tokens[4]
                if len(entry_tokens)>5:
                    x1feat_downstream_of_edit_site_3prime_cp=entry_tokens[-1]
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
                        x1feat_downstream_of_edit_site_3prime_cp=None
                break
    if(to_find!=None):
        for entry in annotation[7::]:
            if entry.startswith(to_find):
                x1feat_downstream_of_edit_site_3prime_cp=entry.strip().split(' ')[-1]
                pos_info=entry.strip().split(' ')[1]
                pos_start=pos_info.split('..')[0]
                pos_end=pos_info.split('..')[1]
                x1feat_downstream_of_edit_site_length_complementary_strand=int(pos_end)-int(pos_start)+1 
                break
    return x1feat_downstream_of_edit_site \
        ,x1feat_downstream_of_edit_site_length_editing_strand \
        ,x1feat_downstream_of_edit_site_length_complementary_strand \
        ,x1feat_downstream_fo_edit_site_5prime_cp \
        ,x1feat_downstream_of_edit_site_3prime_cp \
        ,feat_end

def annotate_structure(editing_levels,bprna_data,approach):
    '''
    For each mutation: 
        mfeat
        mfeat_prev
        mfeat_next
        mfeat_same_as_edit    

    editing_feature
    hairping_length
    stem_length

    x1feat_downstream_of_edit_site,
    x1feat_downstream_of_edit_site_length_editing_strand,
    x1feat_downstream_of_edit_site_length_complementary_strand,
    x1feat_downstream_fo_edit_site_5prime_cp,
    x1feat_downstream_of_edit_site_3prime_cp,
    x1feat_downstream_of_edit_site_endpos

    x2feat_downstream_of_edit_site
    x2feat_downstream_of_edit_site_length_editing_strand
    x2feat_downstream_of_edit_site_length_complementary_strand
    x2feat_downstream_of_edit_site_5prime_cp
    x2feat_downstream_of_edit_site_3prime_cp
    '''
    structure_features=dict() 
    for cur_id in bprna_data:
        annotation=bprna_data[cur_id][approach].split('\n')
        features=annotation[5]
        cur_sequence=annotation[3]
        cur_structure=annotation[4]

        
        #Annotate Editing Site Features
        #IMPORTANT -- DATA IS PROVIDED 1-INDEXED. SHIFT TO 0-INDEX FOR USE IN INDEXING FEATURE ARRAY
        try:
            editing_site=int(editing_levels[cur_id]['site'])-1
        except:
            editing_site=None
        editing_feature=features[editing_site]
        structure_features[cur_id]=dict()
        
        #get the stem length
        stem_length=get_structure_length(features,editing_site,'S') 
        #get the hairpin length
        hairpin_length=get_structure_length(features,editing_site,'H')
        for i in range(len(editing_levels[cur_id]['mut'].keys())):
            mfeat=None
            mfeat_prev=None
            mfeat_next=None
            mfeat_same_as_edit=None
        
            #Annotate current mutation 
            #SHIFT mp to 0-indexing as well
            mp=editing_levels[cur_id]['mut'][i]['mp']
            if mp!=None:
                mp=int(mp)-1
                mfeat,mfeat_prev,mfeat_next,mfeat_same_as_edit=get_bprna_feature_labels(mp,features,editing_site)
            editing_levels[cur_id]['mut'][i]['mfeat']=mfeat
            editing_levels[cur_id]['mut'][i]['mfeat_prev']=mfeat_prev
            editing_levels[cur_id]['mut'][i]['mfeat_next']=mfeat_next
            editing_levels[cur_id]['mut'][i]['mfeat_same_as_edit']=mfeat_same_as_edit
                        
            
        #length of internal loop downstream of editing site, type of loop, closing pair 
        x1feat_downstream_of_edit_site \
            ,x1feat_downstream_of_edit_site_length_editing_strand \
            ,x1feat_downstream_of_edit_site_length_complementary_strand \
            ,x1feat_downstream_fo_edit_site_5prime_cp \
            ,x1feat_downstream_of_edit_site_3prime_cp \
            ,x1feat_downstream_of_edit_site_endpos=get_downstream_nonstem_info(features,editing_site,annotation)
        
        x2feat_downstream_of_edit_site \
            ,x2feat_downstream_of_edit_site_length_editing_strand \
            ,x2feat_downstream_of_edit_site_length_complementary_strand \
            ,x2feat_downstream_of_edit_site_5prime_cp \
            ,x2feat_downstream_of_edit_site_3prime_cp \
            ,x2feat_downstream_of_edit_site_endpos=get_downstream_nonstem_info(features,x1feat_downstream_of_edit_site_endpos,annotation)

        #store all to dict
        structure_features[cur_id]['editing_feature']=editing_feature
        structure_features[cur_id]['stem_length']=stem_length
        structure_features[cur_id]['hairpin_length']=hairpin_length
        structure_features[cur_id]['x1feat_downstream_of_edit_site']=x1feat_downstream_of_edit_site
        structure_features[cur_id]['x1feat_downstream_of_edit_site_length_editing_strand']=x1feat_downstream_of_edit_site_length_editing_strand
        structure_features[cur_id]['x1feat_downstream_of_edit_site_length_complementary_strand']=x1feat_downstream_of_edit_site_length_complementary_strand
        structure_features[cur_id]['x1feat_downstream_of_edit_site_5prime_cp']=x1feat_downstream_fo_edit_site_5prime_cp
        structure_features[cur_id]['x1feat_downstream_of_edit_site_3prime_cp']=x1feat_downstream_of_edit_site_3prime_cp
        structure_features[cur_id]['x2feat_downstream_of_edit_site']=x2feat_downstream_of_edit_site
        structure_features[cur_id]['x2feat_downstream_of_edit_site_length_editing_strand']=x2feat_downstream_of_edit_site_length_editing_strand
        structure_features[cur_id]['x2feat_downstream_of_edit_site_length_complementary_strand']=x2feat_downstream_of_edit_site_length_complementary_strand
        structure_features[cur_id]['x2feat_downstream_of_edit_site_5prime_cp']=x2feat_downstream_of_edit_site_5prime_cp
        structure_features[cur_id]['x2feat_downstream_of_edit_site_3prime_cp']=x2feat_downstream_of_edit_site_3prime_cp
    return structure_features,editing_levels 
    
def write_feature_matrix(editing_levels,structure_dict,outf,source):
    outf=open(outf,'w')
    header=None
    for cur_id in editing_levels.keys():
        mut_info_keys=list(editing_levels[cur_id]['mut'][0].keys())
        struct_info_keys=list(structure_dict[cur_id].keys())
        editing_level=editing_levels[cur_id]['level']
        if header==None:
            header=['source_cur_id','editing_level','num_mutations']+mut_info_keys+struct_info_keys
            header='\t'.join([str(i) for i in header])
            outf.write(header+'\n')
        num_mutations=len(editing_levels[cur_id]['mut'].keys())
        for i in range(num_mutations): 
            mut_info=editing_levels[cur_id]['mut'][i]
            struct_info=structure_dict[cur_id]
            if mut_info['mtype']=="wt":
                num_mutations=0
            outf.write(source+'_'+str(cur_id)+'\t'+str(editing_level)+'\t'+str(num_mutations))
            for keyname in mut_info_keys:
                outf.write('\t'+str(mut_info[keyname]))
            for keyname in struct_info_keys:
                outf.write('\t'+str(struct_info[keyname]))
            outf.write('\n')
                           
def get_editing_info(editing_levels_file,approach):
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
    mp=None
    mref=None
    malt=None
    mtype=None
    adist=None

    if mut.lower().startswith('indel'):
        mtype='indel'
        mp=mut[5::].split('-')[0]
    elif mut.lower().startswith('wt'): 
        #This is wild type
        mtype="wt"
    else:
        mp=mut[:-4]
        adist=editing_site-int(mp)
        mref=mut[-4]
        malt=mut[-1]
        mtype="mismatch"
    return mp,adist,mref,malt,mtype

def get_mut_info(editing_levels_dict):
    '''
    generates the features: 
    mp -- position of mutation 
    adist -- distance from mutation to editing site (1d)
    mref -- reference allele 
    malt -- alternate allele 
    mtype -- mismatch, indel, wt 
    '''
    for cur_id in editing_levels_dict:
        mut_syntax=editing_levels_dict[cur_id]['mutation_syntax'].split(',')
        editing_site=editing_levels_dict[cur_id]['site']        
        #keep track of mutation position information
        mut_dict=dict()
        for i in range(len(mut_syntax)):
            mut_annotations=annotate_mutation(mut_syntax[i],editing_site)
            mut_dict[i]=dict()
            mut_dict[i]['mp']=mut_annotations[0]
            mut_dict[i]['adist']=mut_annotations[1]
            mut_dict[i]['mref']=mut_annotations[2]
            mut_dict[i]['malt']=mut_annotations[3]
            mut_dict[i]['mtype']=mut_annotations[4]
            
        editing_levels_dict[cur_id]['mut']=mut_dict 
    return editing_levels_dict

def annotation_base_in_stem_freq(editing_levels_dict,bprna_data,outf,feat_type):
    entries=editing_levels_dict.keys()
    stem_freq=dict()
    for entry in entries:
        bootstraps=bprna_data[entry]['bootstraps']        
        struct_len=len(bprna_data[entry]['annotation'].split('\n')[5])
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
    
    editing_levels_dict=get_editing_info(editing_levels_file,args.approach)
    
    editing_levels_dict=get_mut_info(editing_levels_dict)
    
    #what fraction of bases in bootstrapped samples are in specified_feature?
    if (args.annotate_bootstraps==True):
        annotation_base_in_stem_freq(editing_levels_dict,bprna_data,outf,'S')
        annotation_base_in_stem_freq(editing_levels_dict,bprna_data,outf,'I')
        annotation_base_in_stem_freq(editing_levels_dict,bprna_data,outf,'B')
        annotation_base_in_stem_freq(editing_levels_dict,bprna_data,outf,'H')
    

    #annotate computational/experime
    structure_dict,editing_levels_dict=annotate_structure(editing_levels_dict,bprna_data,args.approach)
    source=args.source
    if source.__contains__('_'):
        source=source.replace('_','.') 
    write_feature_matrix(editing_levels_dict,structure_dict,outf,source)
    
if __name__=="__main__":
    main()
