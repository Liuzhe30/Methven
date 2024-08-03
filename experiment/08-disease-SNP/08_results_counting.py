# counting results
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
import warnings
warnings.filterwarnings('ignore')

output_path = 'data/pred_merged/'
snp_list = ['rs947474','rs3806624', 'rs7731626', 'rs2234067','rs2233424','rs3824660','rs968567','rs3218251']

model = 'small'
merged_df = pd.DataFrame()
for snp in snp_list:
    data = pd.read_csv(output_path + snp + '_' + model + '.csv')

    # |Δslope|>0.5,slope(t0)>=0
    delta1 = len(data[(data['delta']>0.5)&(data['t0']>0)].reset_index(drop=True))
    
    # |Δslope|<=0.5
    delta2 = len(data[(data['delta']>=-0.5)&(data['delta']<=0.5)].reset_index(drop=True))
    
    # |Δslope|>0.5,slope(t0)<0 
    delta3 = len(data[(data['delta']<-0.5)&(data['t0']<0)].reset_index(drop=True))
    
    # others
    delta4 = len(data) - delta1 - delta2 - delta3

    merged_df = merged_df._append({'SNP':snp,'type':'Up-enhanced','count of CpG sites':delta1},ignore_index=True)
    merged_df = merged_df._append({'SNP':snp,'type':'Unaffected','count of CpG sites':delta2},ignore_index=True)
    merged_df = merged_df._append({'SNP':snp,'type':'Reduced','count of CpG sites':delta4},ignore_index=True)
    merged_df = merged_df._append({'SNP':snp,'type':'Down-enhanced','count of CpG sites':delta3},ignore_index=True)

print(merged_df)
merged_df.to_csv(output_path + 'merged_' + model + '.csv')


model = 'large'
merged_df = pd.DataFrame()
for snp in snp_list:
    data = pd.read_csv(output_path + snp + '_' + model + '.csv')

    # |Δslope|>0.5,slope(t0)>=0
    delta1 = len(data[(data['delta']>0.5)&(data['t0']>0)].reset_index(drop=True))
    
    # |Δslope|<=0.5
    delta2 = len(data[(data['delta']>=-0.5)&(data['delta']<=0.5)].reset_index(drop=True))
    
    # |Δslope|>0.5,slope(t0)<0 
    delta3 = len(data[(data['delta']<-0.5)&(data['t0']<0)].reset_index(drop=True))
    
    # others
    delta4 = len(data) - delta1 - delta2 - delta3

    merged_df = merged_df._append({'SNP':snp,'type':'Up-enhanced','count of CpG sites':delta1},ignore_index=True)
    merged_df = merged_df._append({'SNP':snp,'type':'Unaffected','count of CpG sites':delta2},ignore_index=True)
    merged_df = merged_df._append({'SNP':snp,'type':'Reduced','count of CpG sites':delta4},ignore_index=True)
    merged_df = merged_df._append({'SNP':snp,'type':'Down-enhanced','count of CpG sites':delta3},ignore_index=True)

print(merged_df)
merged_df.to_csv(output_path + 'merged_' + model + '.csv')
