import argparse
import pandas as pd
from statistics import variance 
import pdb

def parse_args():
    parser=argparse.ArgumentParser(description="Rank feature importance across substrates")
    parser.add_argument("--importance_matrix",nargs="+")
    parser.add_argument("--substrates",nargs="+")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    feat_val_dict=dict()
    feat_rank_dict=dict()
    features=set([])
    substrates=list(args.substrates) 
    for i in range(len(args.importance_matrix)):
        cur_substrate=substrates[i]
        f=args.importance_matrix[i] 
        data=pd.read_table(f,sep='\t',header=0)
        data=data.sort_values(by=['%IncMSE'],ascending=False)
        cur_rank=0
        for index,row in data.iterrows():
            feature=index
            features.add(feature) 
            cur_rank+=1
            mse_increase=row['%IncMSE']
            if feature not in feat_val_dict:
                feat_val_dict[feature]=dict()
            if feature not in feat_rank_dict:
                feat_rank_dict[feature]=dict()
            feat_val_dict[feature][cur_substrate]=mse_increase
            feat_rank_dict[feature][cur_substrate]=cur_rank
    print("generated rank dictionaries")
    outf=open(args.outf,'w')
    outf.write('Feature') 
    for substrate in substrates:
        outf.write('\t'+substrate+'_Rank'+'\t'+substrate+"_%IncMSE")
    outf.write('\tMean_Rank\tMean_%IncMSE\tVar_Rank\tVar_%IncMSE\n')
    for feature in features:
        outf.write(feature)
        for substrate in substrates:
            try:
                cur_rank=feat_rank_dict[feature][substrate]
            except:
                cur_rank="NA"
            try:
                cur_value=round(feat_val_dict[feature][substrate],3)
            except:
                cur_value="NA"
            outf.write('\t'+str(cur_rank)+'\t'+str(cur_value))
        all_values=[float(i) for i in feat_val_dict[feature].values()]
        all_ranks=[float(i) for i in feat_rank_dict[feature].values()]
        mean_value=round(sum(all_values)/len(all_values),3)
        mean_rank=round(sum(all_ranks)/len(all_ranks),3)
        if len(all_values)>1:
            var_value=round(variance(all_values),3) 
            var_rank=round(variance(all_ranks),3)
        else:
            var_value=0
            var_rank=0 
        outf.write('\t'+'\t'.join([str(mean_rank),str(mean_value),str(var_rank),str(var_value)])+'\n')
        
if __name__=="__main__":
    main()
    
    
