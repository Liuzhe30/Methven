# merge dnabert embedding & downsampling large-datasets

import pandas as pd
import numpy as np
from sklearn.utils import shuffle
pd.set_option('display.max_columns', None)
import warnings
warnings.filterwarnings('ignore')

atac_path = '../../datasets/atac_cut/'
atac_ori_path = '../../datasets/atac_mapping_post/'
dnabert_path_ori = '../../datasets/dnabert_embedding/'
dnabert_path = '../../datasets/dnabert_embedding_large/'
merge_path = '../../datasets/merge_dnabert/'
output_path1 = '../../datasets/final/dnabert/'
output_path2 = '../../datasets/final/origin/'
clean_path = '../../data/original_split/'

model = 'large'

# down-sampling datasets
for i in range(22):
    chr = 'chr' + str(i+1)
    data = pd.read_pickle(dnabert_path_ori + chr + '_' + model + '.dataset')
    data = shuffle(data)
    sample_num = data.shape[0]
    new_data = data.sample(n=int(sample_num/2),replace=False).reset_index(drop=True)
    new_data.to_pickle(dnabert_path + chr + '_' + model + '.dataset')

# merge embedding & atac-seq
for i in range(22):
    chr = 'chr' + str(i+1)
    data_atac = pd.read_pickle(atac_path + chr + '_' + model + '.dataset')
    data_atac_ori = pd.read_pickle(atac_ori_path + chr + '_' + model + '.dataset')
    data_dnabert = pd.read_pickle(dnabert_path + chr + '_' + model + '.dataset')
    data_dnabert['atac_between'] = 0
    data_dnabert['atac_between'] = data_dnabert['atac_between'].astype('object')
    data_dnabert['atac_between_ori'] = 0
    data_dnabert['atac_between_ori'] = data_dnabert['atac_between_ori'].astype('object')
    for i in range(len(data_dnabert)):
        cpg = data_dnabert['CpG'][i]
        snp = data_dnabert['SNP'][i]
        data_dnabert['atac_between'][i] = data_atac[(data_dnabert['CpG']==cpg)&(data_atac['SNP']==snp)]['atac_between']
        data_dnabert['atac_between_ori'][i] = data_atac_ori[(data_atac_ori['CpG']==cpg)&(data_atac_ori['SNP']==snp)]['atac_between']
    data_dnabert.to_pickle(merge_path + chr + '_' + model + '.dataset')

# merge datasets
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

train_dnabert = train_data[['CpG','SNP','Beta','Ref','Alt','CHR','CpG_POS','SNP_POS','label','distance','seq_len','atac_between','dnabert_before','dnabert_after']]
valid_dnabert = valid_data[['CpG','SNP','Beta','Ref','Alt','CHR','CpG_POS','SNP_POS','label','distance','seq_len','atac_between','dnabert_before','dnabert_after']]
test_dnabert = test_data[['CpG','SNP','Beta','Ref','Alt','CHR','CpG_POS','SNP_POS','label','distance','seq_len','atac_between','dnabert_before','dnabert_after']]
train_dnabert.to_pickle(output_path1 + model + '_train.dataset')
valid_dnabert.to_pickle(output_path1 + model + '_valid.dataset')
test_dnabert.to_pickle(output_path1 + model + '_test.dataset')

train_origin = train_data[['CpG','SNP','Beta','Ref','Alt','CHR','CpG_POS','SNP_POS','label','distance','seq_len','atac_between_ori','seq_before','seq_after']]
valid_origin = valid_data[['CpG','SNP','Beta','Ref','Alt','CHR','CpG_POS','SNP_POS','label','distance','seq_len','atac_between_ori','seq_before','seq_after']]
test_origin = test_data[['CpG','SNP','Beta','Ref','Alt','CHR','CpG_POS','SNP_POS','label','distance','seq_len','atac_between_ori','seq_before','seq_after']]
train_origin.to_pickle(output_path2 + model + '_train.dataset')
valid_origin.to_pickle(output_path2 + model + '_valid.dataset')
test_origin.to_pickle(output_path2 + model + '_test.dataset')

train_clean = train_data[['CpG','SNP','Beta','Ref','Alt','CHR','CpG_POS','SNP_POS']]
valid_clean = valid_data[['CpG','SNP','Beta','Ref','Alt','CHR','CpG_POS','SNP_POS']]
test_clean = test_data[['CpG','SNP','Beta','Ref','Alt','CHR','CpG_POS','SNP_POS']]
train_clean.to_pickle(clean_path + model + '_train.dataset')
valid_clean.to_pickle(clean_path + model + '_valid.dataset')
test_clean.to_pickle(clean_path + model + '_test.dataset')