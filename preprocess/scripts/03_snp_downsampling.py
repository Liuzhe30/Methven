import os
import pandas as pd
import numpy as np
import subprocess
from tqdm import tqdm
import random
import json

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

#'''
# step 1: split by chr
for i in range(22):
    sub_table = cross_table[cross_table['CHR'] == str(i+1)].reset_index(drop = True)
    sub_table.to_pickle('../../datasets/middlefile/meQTL_annotation_CpG_SNP_chr/chr' + str(i+1) + '.pkl')
#'''

#'''
# step 2: downsample SNPs
down_snp_dict = {}
for i in range(22):
    sub_table = pd.read_pickle('../../datasets/middlefile/meQTL_annotation_CpG_SNP_chr/chr' + str(i+1) + '.pkl')
    snp_list = sub_table['SNP'].unique().tolist()
    #print(len(snp_list))
    samples = random.sample(snp_list, 2000)
    down_snp_dict['chr' + str(i+1)] = samples
    json_str = json.dumps(down_snp_dict, indent=4)
with open('../../datasets/middlefile/downsample_SNP.json', 'w') as json_file:
    json_file.write(json_str)
#'''   
