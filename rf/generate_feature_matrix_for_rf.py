import pandas as pd
import pdb
editing_levels=pd.read_csv("Neil1_Editing_Level_For_Anna_2018-01-23.csv",header=0,sep='\t')
bprna_output_dir="/srv/scratch/annashch/adar_editing/data_from_xin/Neil1_raw_data-bpRNA-results"
feature_dict=dict()
start_struct=1
end_struct=293
editing_site=50
outf=open("rf_features.txt",'w')
outf.write('\t'.join(["ave_editing_level",
                      "edit_feat",
                      "prev_feat",
                      "next_feat",
                      "mp1",
                      "mref1",
                      "malt1",
                      "mtype1",
                      "mfeat1",
                      "mfeat1_prev",
                      "mfeat1_next",
                      "adist1",
                      "mp2",
                      "mref2",
                      "malt2",
                      "mtype2",
                      "mfeat2",
                      "mfeat2_prev",
                      "mfeat2_next",
                      "adist2"])+'\n')

for index,row in editing_levels.iterrows():
    seq_index=row['RNA_ID']
    if seq_index<10:
        seq_index="00"+str(seq_index)
    elif seq_index<100:
        seq_index="0"+str(seq_index)
    else:
        seq_index=str(seq_index) 
    bprna_data=open(bprna_output_dir+'/'+str(seq_index)+"-reactivity.dot.st",'r').read().strip().split('\n')
    features=bprna_data[5]
    edit_feat=features[editing_site-1]
    #previous feature
    for i in range(editing_site-2,0,-1):
        if features[i]!=edit_feat:
            prev_feat=features[i]
            break
    #subsequent feature
    for i in range(editing_site,len(features)):
        if features[i]!=edit_feat:
            next_feat=features[i]
            break
    ave_editing_level=row["Ave_Editing_Level"]
    mut_syntax=row["Mut_Syntax"].split(',')
    if len(mut_syntax) not in [1,2]:
        pdb.set_trace()
    mp1=None
    mp2=None
    mref1=None
    mref2=None
    malt1=None
    malt2=None
    mtype1=None
    mtype2=None
    mfeat1=None
    mfeat1_prev=None
    mfeat1_next=None
    mfeat2=None
    mfeat2_prev=None
    mfeat2_next=None
    adist1=None
    adist2=None
    if len(mut_syntax)==2:
        #annotate the first mutation
        first_mut=mut_syntax[0]
        mp1=first_mut[:-4]
        adist1=editing_site-int(mp1)
        mref1=first_mut[-4]
        malt1=first_mut[-1]
        mtype1="mismatch"
        mfeat1=features[int(mp1)-1]

        for i in range(int(mp1)-2,0,-1):
            if features[i]!=mfeat1:
                mfeat1_prev=features[i]
                break
            #subsequent feature
        for i in range(int(mp1),len(features)):
            if features[i]!=mfeat1:
                mfeat1_next=features[i]
                break
        
        second_mut=mut_syntax[1]
        mp2=second_mut[:-4]
        adist2=editing_site-int(mp2)
        mref2=second_mut[-4]
        malt2=second_mut[-1]
        mtype2="mismatch"
        mfeat2=features[int(mp2)-1]

        for i in range(int(mp2)-2,0,-1):
            if features[i]!=mfeat2:
                mfeat2_prev=features[i]
                break
            #subsequent feature
        for i in range(int(mp2),len(features)):
            if features[i]!=mfeat1:
                mfeat1_next=features[i]
                break

        
    elif "to" in mut_syntax[0]:
        first_mut=mut_syntax[0]
        mp1=first_mut[:-4]
        adist1=editing_site-int(mp1)
        mref1=first_mut[-4]
        malt1=first_mut[-1]
        mtype1="mismatch"
        mfeat1=features[int(mp1)-1]
        for i in range(int(mp1)-2,0,-1):
            if features[i]!=mfeat1:
                mfeat1_prev=features[i]
                break
            #subsequent feature
        for i in range(int(mp1),len(features)):
            if features[i]!=mfeat1:
                mfeat1_next=features[i]
                break        
    elif "Indel" in mut_syntax[0]:
        mtype1="indel"
    elif "WT" in mut_syntax[0]:
        mtype1="wt"
    else:
        pdb.set_trace() 
    #write line to output file
    outf.write('\t'.join([str(i) for i in [ave_editing_level,
                          edit_feat,
                          prev_feat,
                          next_feat,
                          mp1,
                          mref1,
                          malt1,
                          mtype1,
                          mfeat1,
                          mfeat1_prev,
                          mfeat1_next,
                          adist1,
                          mp2,
                          mref2,
                          malt2,
                          mtype2,
                          mfeat2,
                          mfeat2_prev,
                          mfeat2_next,
    adist2]])+'\n')
