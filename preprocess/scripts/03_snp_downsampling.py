import os
import pandas as pd
import numpy as np
import subprocess
from tqdm import tqdm
import random
import json
import warnings
warnings.filterwarnings('ignore')

cross_path = '../../datasets/middlefile/meQTL_annotation_CpG_SNP.pkl'
meQTL_path = '../../datasets/middlefile/single_meQTL.txt'

cross_table = pd.read_pickle(cross_path)
#print(cross_table.head())
'''
          CpG          SNP   Beta Ref Alt CHR   CpG_POS   SNP_POS
0  cg11913416    rs1262461 -1.179   A   G   1  13582516  13581610
1  cg25722041  rs114812780  1.233   G   T   1   8563413   8468794
2  cg25722041   rs12568293  1.233   T   C   1   8563413   8480852
3  cg25722041   rs12117910  1.233   A   C   1   8563413   8496893
4  cg25722041   rs12403339  1.233   G   A   1   8563413   8498232
'''

'''
# step 1: calculate distance
cross_table = cross_table.reset_index(drop=True)
cross_table['distance'] = 0
for i in tqdm(range(cross_table.shape[0])):
    cpg_pos = cross_table['CpG_POS'][i]
    snp_pos = cross_table['SNP_POS'][i]
    distance = np.abs(int(cpg_pos) - int(snp_pos))
    cross_table['distance'][i] = distance
cross_table.to_pickle('../../datasets/middlefile/meQTL_annotation_CpG_SNP_distance.pkl')
print(cross_table.head())    
'''
'''
          CpG          SNP   Beta Ref Alt CHR   CpG_POS   SNP_POS  distance
0  cg11913416    rs1262461 -1.179   A   G   1  13582516  13581610       906
1  cg25722041  rs114812780  1.233   G   T   1   8563413   8468794     94619
2  cg25722041   rs12568293  1.233   T   C   1   8563413   8480852     82561
3  cg25722041   rs12117910  1.233   A   C   1   8563413   8496893     66520
4  cg25722041   rs12403339  1.233   G   A   1   8563413   8498232     65181
'''

'''
# step 2: split by chr
for i in range(22):
    sub_table = cross_table[cross_table['CHR'] == str(i+1)].reset_index(drop = True)
    sub_table.to_pickle('../../datasets/middlefile/meQTL_annotation_CpG_SNP_chr/chr' + str(i+1) + '.pkl')
'''

# step 3: downsample SNPs
# orignal ratio: small:middle:large = 1:5:10
down_snp_dict = {}
for i in range(22):
    down_snp_dict['chr' + str(i+1)] = {}
    sub_table = pd.read_pickle('../../datasets/middlefile/meQTL_annotation_CpG_SNP_chr/chr' + str(i+1) + '.pkl')

    len_table = sub_table[sub_table['distance'] <= 1_000]
    snp_list = len_table['SNP'].unique().tolist()
    print(len(snp_list))
    samples = random.sample(snp_list, 2000)
    down_snp_dict['chr' + str(i+1)]['small'] = samples

    len_table = sub_table[(sub_table['distance'] > 1_000) & (sub_table['distance'] <= 10_000)]
    snp_list = len_table['SNP'].unique().tolist()
    print(len(snp_list))
    samples = random.sample(snp_list, 2000)
    down_snp_dict['chr' + str(i+1)]['middle'] = samples

    len_table = sub_table[(sub_table['distance'] > 10_000) & (sub_table['distance'] <= 100_000)]
    snp_list = len_table['SNP'].unique().tolist()
    print(len(snp_list))
    samples = random.sample(snp_list, 2000)
    down_snp_dict['chr' + str(i+1)]['large'] = samples

    len_table = sub_table[(sub_table['distance'] > 100_000) & (sub_table['distance'] <= 1_000_000)]
    snp_list = len_table['SNP'].unique().tolist()
    print(len(snp_list))
    samples = random.sample(snp_list, 2000)
    down_snp_dict['chr' + str(i+1)]['huge'] = samples


json_str = json.dumps(down_snp_dict, indent=4)
with open('../../datasets/middlefile/downsample_SNP.json', 'w') as json_file:
    json_file.write(json_str)
  
