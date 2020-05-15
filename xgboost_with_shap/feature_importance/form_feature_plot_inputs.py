import pdb 
import pandas as pd
import numpy as np
neil1=pd.read_csv("NEIL1_comp_shap_mean_abs.txt",header=0,sep='\t')
ttyh2=pd.read_csv("TTYH2_comp_shap_mean_abs.txt",header=0,sep='\t')
ajuba=pd.read_csv("AJUBA_comp_shap_mean_abs.txt",header=0,sep='\t')
neil1['substrate']="NEIL1"
ttyh2['substrate']="TTYH2"
ajuba['substrate']="AJUBA" 
data=pd.concat([neil1,ttyh2,ajuba])
data.to_csv("FeatureSHAP.tsv",header=True,index=False,sep='\t')

data_dict=dict()
for index,row in data.iterrows():
    feature=row['feature']
    meanshap=row['mean_abs_shap']
    substrate=row['substrate'] 
    if feature not in data_dict:
        data_dict[feature]=dict()
    data_dict[feature][substrate]=float(meanshap)
data=pd.DataFrame(data_dict)
data=data.transpose()
#drop any features not in all 3 substrates
data=data.dropna()
print(data.head())
data=data/np.sum(data,axis=0)
data=data*100
print(data.head())
mean_impact=data.mean(axis=1)
std_impact=data.std(axis=1)
data['MeanImpact']=mean_impact
data['StdImpact']=std_impact
data['Rank']=data['MeanImpact'].rank(ascending=False)
data.to_csv('FeatureRanksByShap.tsv',sep='\t',index=True,header=True)



