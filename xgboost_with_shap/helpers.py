#helper functions for training xgboost models to predict Adar editing levels
import random 
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder,OneHotEncoder
from scipy.stats import spearmanr,pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error, average_precision_score, roc_auc_score
import math 


def mape(labels,preds,pseudocount=0.01):
    return np.mean(np.abs((labels - preds) /(pseudocount+labels)))

def print_aggregate_peformance_metrics(train_labels,train_preds,eval_labels,eval_preds,test_labels,test_preds,subset=None,round_places=3):
    metrics=aggregate_performance_metrics(train_labels,train_preds,eval_labels,eval_preds,test_labels,test_preds,subset=subset)
    for key in metrics:
        print(key)
        for subkey in metrics[key]:
            cur_metric=metrics[key][subkey]
            print('\t'+subkey+' : '+str(cur_metric))
            

def aggregate_performance_metrics(train_labels,train_preds,eval_labels,eval_preds,test_labels,test_preds,subset=None):
    if subset!=None:
        train_labels=train_labels[subset]
        tain_preds=train_preds[subset]
        eval_labels=eval_labels[subset]
        eval_preds=eval_preds[subset]
        test_labels=test_labels[subset]
        test_preds=test_preds[subset]
    metrics=dict()     
    metrics['Train']=get_performance_metrics(train_labels,train_preds)
    metrics['Eval']=get_performance_metrics(eval_labels,eval_preds)
    metrics['Test']=get_performance_metrics(test_labels,test_preds)
    return metrics 


def get_performance_metrics(labels,preds):
    metrics=dict() 
    metrics['Spearman corr']=spearmanr(labels,preds)
    metrics['Pearson corr']=pearsonr(labels,preds)
    metrics['MAE']=mean_absolute_error(labels,preds)
    metrics['MAPE']=mape(labels,preds)
    metrics['RMSE']=math.sqrt(mean_squared_error(labels,preds))
    metrics['auPRC']=average_precision_score(labels>=0.5,preds)
    metrics['auROC']=roc_auc_score(labels>=0.5,preds)
    return metrics



def shap_contribs_subgroup(shap_values,feature_group):
    shap_mean_abs_by_subject=abs(shap_values).mean(axis=0)
    total_shap=shap_mean_abs_by_subject.sum()
    shap_mean_abs_by_subject_norm=shap_mean_abs_by_subject/total_shap
    shap_subset=shap_mean_abs_by_subject_norm[filter_features(shap_values,feature_group)]
    return shap_subset.sum()

    

def filter_features(df,feature_group):
    filtered_columns=[] 
    for f in feature_group:
        new_cols=list(df.filter(regex=f,axis=1).columns)
        filtered_columns=filtered_columns+new_cols
    return filtered_columns

#Identify any columns that contain features with constant value across datapoints
def get_singleval_features(X):
    todrop=[]
    for column in X.columns:
        if len(set(X[column]))==1:
            todrop.append(column)
    return todrop 
#Identify any columns that contain features with null values across datapoints 
def get_all_null_features(X):
    all_null=[]
    for c in X.columns:
        if (pd.isna(X[c]).all())==True:
            all_null.append(c)
    return all_null

def format_for_xgboost(X):
    first=True
    for i in range(X.shape[1]):
        feature=X[X.columns[i]]
        feat_name=X.columns[i]
        if feature.dtype not in [float,int]:
            #We need to handle NA values 
            feature[pd.isna(feature)]="NA"
            label_encoder = LabelEncoder()
            feature = label_encoder.fit_transform(feature)
            feature = feature.reshape(X.shape[0], 1)
            onehot_encoder = OneHotEncoder(sparse=False)
            feature = onehot_encoder.fit_transform(feature)
            #Drop the column corresponding to NA
            feature=pd.DataFrame(feature,columns=[feat_name+":"+str(j) for j in label_encoder.classes_])
            embedding_map = dict(zip(label_encoder.classes_, label_encoder.transform(label_encoder.classes_)))
            if "NA" in embedding_map:
                todrop=feat_name+":NA"
                feature=feature.drop([todrop],axis=1)
        else: 
            feature = feature.values.reshape(X.shape[0], 1)
            feature=pd.DataFrame(feature,columns=[feat_name])
        if first==True: 
            transformed=feature
            first=False
        else: 
            transformed=pd.concat((transformed,feature),axis=1)
    print(transformed.shape)
    return transformed

def split_train_test_eval_by_mut_pos(data,train_split_percent=0.70,eval_split_percent=0.15,test_split_percent=0.15):
    valcount,positions=np.histogram(data['mut_pos'], bins=np.unique(np.sort(data['mut_pos'])))
    pos_to_count= dict(zip(positions,valcount))
    print(pos_to_count)
    positions=list(pos_to_count.keys())
    total=data.shape[0] 
    random.shuffle(positions)
    added=0
    train_pos=[]
    eval_pos=[] 
    test_pos=[] 
    num_train=total*train_split_percent
    num_eval=total*eval_split_percent 
    num_test=total*test_split_percent 
    for pos in positions: 
        added+=pos_to_count[pos]
        if added < num_train: 
            train_pos.append(pos)
        elif added < (num_train+num_eval): 
            eval_pos.append(pos) 
        else: 
            test_pos.append(pos)
    train_split=data[data['mut_pos'].isin(train_pos)]
    eval_split=data[data['mut_pos'].isin(eval_pos)]
    test_split=data[data['mut_pos'].isin(test_pos)]
    return train_split,eval_split,test_split
