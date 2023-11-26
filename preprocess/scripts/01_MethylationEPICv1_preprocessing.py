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
print(new_data.head())

'''
       IlmnID CHR_hg38   Start_hg38     End_hg38
0  cg07881041    chr19    5236004.0    5236006.0
1  cg18478105    chr20   63216297.0   63216299.0
2  cg23229610     chr1    6781064.0    6781066.0
3  cg03513874     chr2  197438741.0  197438743.0
'''

# save data
new_data.to_pickle(save_path + 'cpg_all.pkl')

chr_list = new_data['CHR_hg38'].unique()
for i in chr_list:
    chr_data = new_data[new_data['CHR_hg38'] == i]
    chr_data = chr_data.reset_index(drop=True)
    print(chr_data)
    chr_data.to_pickle(save_path + 'cpg_' + i + '.pkl')