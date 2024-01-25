import os
import pandas as pd
import numpy as np

file_path = '../../datasets/raw/infinium-methylationepic-v-1-0-b5-manifest-file.csv' # remove head and [control] items
save_path = '../../datasets/middlefile/clean_epic/'

# load raw data
data = pd.read_csv(file_path)
# print(data.shape) # (865918, 52)

# filter columns
new_data = data[['IlmnID','CHR_hg38','Start_hg38','End_hg38']]
new_data.dropna(axis=0, how='any', inplace=True)
new_data[['Start_hg38']] = new_data[['Start_hg38']].astype(int)
new_data[['End_hg38']] = new_data[['End_hg38']].astype(int)
print(new_data.head())
print(new_data.shape) # (865880, 4)

'''
       IlmnID CHR_hg38  Start_hg38   End_hg38
0  cg07881041    chr19     5236004    5236006
1  cg18478105    chr20    63216297   63216299
2  cg23229610     chr1     6781064    6781066
3  cg03513874     chr2   197438741  197438743
4  cg09835024     chrX    24054522   24054524
'''

# save data
new_data.to_pickle(save_path + 'cpg_all.pkl')

chr_list = new_data['CHR_hg38'].unique()
for i in chr_list:
    chr_data = new_data[new_data['CHR_hg38'] == i]
    chr_data = chr_data.reset_index(drop=True)
    print(chr_data)
    chr_data.to_pickle(save_path + 'cpg_' + i + '.pkl')