import pandas as pd
import numpy as np
neil1=pd.read_csv("neil1_computational.features.csv",header=0,sep=',')
ttyh2_bc=pd.read_csv("ttyh2_bc_computational.features.csv",header=0,sep=',')
ttyh2_ecs=pd.read_csv("ttyh2_ecs_bc_computational.features.csv",header=0,sep=',')
ajuba=pd.read_csv("ajuba_bc_computational.features.csv",header=0,sep=',')
outf=open("mut_and_edit.txt",'w')
outf.write("ID\tMut\tEditing\n")
neil_editing=47
ttyh2_bc_editing=87
ttyh2_ecs_editing=87
ajuba_editing=74
for index,row in neil1.iterrows():
    cur_id="NEIL1_"+str(row['rna_id'])
    cur_mut=row['mut_pos']
    if np.isnan(cur_mut):
        cur_mut=-100
    cur_edit=neil_editing
    outf.write(cur_id+'\t'+str(cur_mut)+'\t'+str(cur_edit)+'\n')
for index,row in ttyh2_bc.iterrows():
    cur_id="TTYH2_BC_"+str(row['rna_id'])
    cur_mut=row['mut_pos']
    if np.isnan(cur_mut):
        cur_mut=-100
    cur_edit=ttyh2_bc_editing
    outf.write(cur_id+'\t'+str(cur_mut)+'\t'+str(cur_edit)+'\n')
for index,row in ttyh2_ecs.iterrows():
    cur_id="TTYH2_ECS_"+str(row['rna_id'])
    cur_mut=row['mut_pos']
    if np.isnan(cur_mut):
        cur_mut=-100
    cur_edit=ttyh2_ecs_editing
    outf.write(cur_id+'\t'+str(cur_mut)+'\t'+str(cur_edit)+'\n')
for index,row in ajuba.iterrows():
    cur_id="AJUBA_"+str(row['rna_id'])
    cur_mut=row['mut_pos']
    if np.isnan(cur_mut):
        cur_mut=-100
    cur_edit=ajuba_editing
    outf.write(cur_id+'\t'+str(cur_mut)+'\t'+str(cur_edit)+'\n')
    
outf.close()


