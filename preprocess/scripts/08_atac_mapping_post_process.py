# post process the 'nan' in ATAC-seq
import numpy as np
import pandas as pd
import math
pd.set_option('display.max_columns', None)

def check(check_list):
    for elem in check_list:
        if(float(elem) == np.nan):
            return 0
        if(math.isnan(elem) == True):
            return 0
    return 1

file_path = '../../datasets/atac_mapping/'
output_path = '../../datasets/atac_mapping_post/'

model_size = {'small':1_000,'middle':10_000,'large':100_000,'huge':1_000_000}
for i in range(22):
    chr = 'chr' + str(i+1)
    for model in model_size.keys():
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        data_check = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        for i in range(len(data_check)):
            atac_between = data_check['atac_between'][i]
            cpg = data_check['CpG'][i]
            snp = data_check['SNP'][i]
            if(check(atac_between)==0):
                data = data.drop(data[(data['CpG']==cpg)&(data['SNP']==snp)].index)
        data = data.reset_index(drop=True)
        data.to_pickle(output_path + chr + '_' + model + '.dataset')