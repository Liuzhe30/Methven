import os
import pandas as pd
import numpy as np

file_path = '../../datasets/raw/top_meQTL.txt'
save_path = '../../datasets/middlefile/clean_meQTL.pkl'

# load raw data
data = pd.read_csv(file_path, sep='\t')
print(data.shape) # (249710, 17)

# 1 remove trans-meQTLs
data = data[~(data['Cis/Trans'] == 'trans')]
print(data.shape) # (244491, 17)

# 2 remove multi-point variants
data = data[~(data['A1'].str.len() > 1)]
data = data[~(data['A2'].str.len() > 1)]
print(data.shape) # (199960, 17)

print(data.head())

'''
          CpG       Top SNP  CpG chr   CpG pos  SNP chr  SNP pos A1 A2    MAF      Beta        SE              P  FDR  N     n Effects Cis/Trans
0  cg06325811  1:753405_C_A        1  796328.0        1   753405  C  A  0.158  0.371869  0.049974   1.040000e-13  0.0  2  1531   ?+?+?       cis
1  cg16619049  1:798959_G_A        1  805541.0        1   798959  G  A  0.214  0.769962  0.061491   6.230000e-36  0.0  2  1545   ?++??       cis
3  cg12445832  1:834999_G_A        1  834295.0        1   834999  G  A  0.203 -0.386778  0.048310   1.240000e-15  0.0  2  1545   ?--??       cis
5  cg08128007  1:838387_T_C        1  839435.0        1   838387  T  C  0.205 -1.044016  0.039710  2.850000e-152  0.0  2  1545   ?--??       cis
6  cg23733394  1:836529_C_G        1  839752.0        1   836529  C  G  0.198 -1.163682  0.040631  2.480000e-180  0.0  2  1545   ?--??       cis
'''

# save data
data.to_pickle(save_path)