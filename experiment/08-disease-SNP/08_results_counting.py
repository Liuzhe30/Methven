# counting results
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
import warnings
warnings.filterwarnings('ignore')

output_path = 'data/pred_merged/'
snp_list = ['rs3806624', 'rs7731626', 'rs2234067','rs2233424','rs947474','rs3824660','rs968567','rs3218251']

model = 'small'
merged_df = pd.DataFrame()
for snp in snp_list:
    data = pd.read_csv(output_path + snp + '_' + model + '.csv')
    delta1 = len(data[data['delta']<-0.5].reset_index(drop=True))
    delta3 = len(data[data['delta']>0.5].reset_index(drop=True))
    delta2 = len(data) - delta1 - delta3
    merged_df = merged_df._append({'SNP':snp,'type':'Δslope<-0.5','count':delta1},ignore_index=True)
    merged_df = merged_df._append({'SNP':snp,'type':'|Δslope|<=0.5','count':delta2},ignore_index=True)
    merged_df = merged_df._append({'SNP':snp,'type':'Δslope>0.5','count':delta3},ignore_index=True)
print(merged_df)
merged_df.to_csv(output_path + 'merged_' + model + '.csv')
'''
          SNP           type  count
0   rs3806624    Δslope<-0.5     11
1   rs3806624  |Δslope|<=0.5     32
2   rs3806624     Δslope>0.5     19
3   rs7731626    Δslope<-0.5      2
4   rs7731626  |Δslope|<=0.5      9
5   rs7731626     Δslope>0.5      1
6   rs2234067    Δslope<-0.5      0
7   rs2234067  |Δslope|<=0.5     26
8   rs2234067     Δslope>0.5      0
9   rs2233424    Δslope<-0.5     18
10  rs2233424  |Δslope|<=0.5     12
11  rs2233424     Δslope>0.5     12
12   rs947474    Δslope<-0.5      2
13   rs947474  |Δslope|<=0.5      2
14   rs947474     Δslope>0.5      3
15  rs3824660    Δslope<-0.5     32
16  rs3824660  |Δslope|<=0.5     17
17  rs3824660     Δslope>0.5     18
18   rs968567    Δslope<-0.5      6
19   rs968567  |Δslope|<=0.5     17
20   rs968567     Δslope>0.5      9
21  rs3218251    Δslope<-0.5      4
22  rs3218251  |Δslope|<=0.5     11
23  rs3218251     Δslope>0.5      5
'''

model = 'large'
merged_df = pd.DataFrame()
for snp in snp_list:
    data = pd.read_csv(output_path + snp + '_' + model + '.csv')
    delta1 = len(data[data['delta']<-0.5].reset_index(drop=True))
    delta3 = len(data[data['delta']>0.5].reset_index(drop=True))
    delta2 = len(data) - delta1 - delta3
    merged_df = merged_df._append({'SNP':snp,'type':'Δslope<-0.5','count':delta1},ignore_index=True)
    merged_df = merged_df._append({'SNP':snp,'type':'|Δslope|<=0.5','count':delta2},ignore_index=True)
    merged_df = merged_df._append({'SNP':snp,'type':'Δslope>0.5','count':delta3},ignore_index=True)
print(merged_df)
merged_df.to_csv(output_path + 'merged_' + model + '.csv')

'''
          SNP           type  count
0   rs3806624    Δslope<-0.5      4
1   rs3806624  |Δslope|<=0.5     36
2   rs3806624     Δslope>0.5      4
3   rs7731626    Δslope<-0.5      6
4   rs7731626  |Δslope|<=0.5     43
5   rs7731626     Δslope>0.5     16
6   rs2234067    Δslope<-0.5      2
7   rs2234067  |Δslope|<=0.5     84
8   rs2234067     Δslope>0.5      1
9   rs2233424    Δslope<-0.5     56
10  rs2233424  |Δslope|<=0.5     76
11  rs2233424     Δslope>0.5     52
12   rs947474    Δslope<-0.5      8
13   rs947474  |Δslope|<=0.5     40
14   rs947474     Δslope>0.5      8
15  rs3824660    Δslope<-0.5     21
16  rs3824660  |Δslope|<=0.5     57
17  rs3824660     Δslope>0.5     18
18   rs968567    Δslope<-0.5     60
19   rs968567  |Δslope|<=0.5    124
20   rs968567     Δslope>0.5     47
21  rs3218251    Δslope<-0.5     51
22  rs3218251  |Δslope|<=0.5     87
23  rs3218251     Δslope>0.5     38
'''