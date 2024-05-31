# merge dnabert embedding
import pandas as pd
import numpy as np
from sklearn.utils import shuffle
pd.set_option('display.max_columns', None)
import warnings
warnings.filterwarnings('ignore')

atac_path = '../../datasets/atac_cut/'
dnabert_path = '../../datasets/dnabert_embedding/'
merge_path = '../../datasets/merge_dnabert/'
output_path = '../../datasets/final/dnabert/'

model_size = {'small':10_000,'large':100_000}

# merge embedding & atac-seq
for i in range(22):
    chr = 'chr' + str(i+1)
    for model in model_size.keys():
        data_atac = pd.read_pickle(atac_path + chr + '_' + model + '.dataset')
        data_dnabert = pd.read_pickle(dnabert_path + chr + '_' + model + '.dataset')
        data_atac['dnabert_before'] = 0
        data_atac['dnabert_before'] = data_atac['dnabert_before'].astype('object')
        data_atac['dnabert_after'] = 0
        data_atac['dnabert_after'] = data_atac['dnabert_after'].astype('object')
        for i in range(len(data_atac)):
            cpg = data_atac['CpG'][i]
            snp = data_atac['SNP'][i]
            data_atac['dnabert_before'][i] = data_dnabert[(data_dnabert['CpG']==cpg)&(data_dnabert['SNP']==snp)]['dnabert_before']
            data_atac['dnabert_after'][i] = data_dnabert[(data_dnabert['CpG']==cpg)&(data_dnabert['SNP']==snp)]['dnabert_after']
        data_atac.to_pickle(merge_path + chr + '_' + model + '.dataset')

# merge datasets
for model in model_size.keys():
    model_merged = pd.DataFrame()
    for i in range(22):
        chr = 'chr' + str(i+1)
        data = pd.read_pickle(merge_path + chr + '_' + model + '.dataset')
        model_merged = pd.concat([model_merged, data])
    model_merged = shuffle(model_merged)
    model_merged = model_merged.reset_index(drop=True)
    train_data = model_merged[0:int(0.8*len(model_merged))].reset_index(drop=True)
    valid_data = model_merged[int(0.8*len(model_merged)):int(0.9*len(model_merged))].reset_index(drop=True)
    test_data = model_merged[int(0.9*len(model_merged)):].reset_index(drop=True)
    train_data.to_pickle(output_path + model + '_train.dataset')
    valid_data.to_pickle(output_path + model + '_valid.dataset')
    test_data.to_pickle(output_path + model + '_test.dataset')