# merge dnabert embedding
import pandas as pd
import numpy as np
from sklearn.utils import shuffle
pd.set_option('display.max_columns', None)
import warnings
warnings.filterwarnings('ignore')

cell_list = ['B_Naive','CD4_memory','CD8_memory']
cell = cell_list[0]

atac_path = '../../datasets/atac_cut_eQTL/' + cell
atac_ori_path = '../../datasets/atac_mapping_post_eQTL/' + cell
dnabert_path = '../../datasets/dnabert_embedding_eQTL/' + cell
merge_path = '../../datasets/merge_dnabert_eQTL/' + cell
output_path1 = '../../datasets/final_eQTL/' + cell

model_size = {'small':10_000,'large':100_000}

for model in model_size.keys():
    # merge embedding & atac-seq
    for i in range(22):
        chr = 'chr' + str(i+1)
        data_atac = pd.read_pickle(atac_path + chr + '_' + model + '.dataset')
        data_atac_ori = pd.read_pickle(atac_ori_path + chr + '_' + model + '.dataset')
        data_dnabert = pd.read_pickle(dnabert_path + chr + '_' + model + '.dataset')
        data_atac['dnabert_before'] = 0
        data_atac['dnabert_before'] = data_atac['dnabert_before'].astype('object')
        data_atac['dnabert_after'] = 0
        data_atac['dnabert_after'] = data_atac['dnabert_after'].astype('object')
        data_atac['atac_between_ori'] = 0
        data_atac['atac_between_ori'] = data_atac['atac_between_ori'].astype('object')
        for i in range(len(data_atac)):
            cpg = data_atac['POS'][i]
            snp = data_atac['TSS'][i]
            data_atac['dnabert_before'][i] = data_dnabert[(data_dnabert['POS']==cpg)&(data_dnabert['TSS']==snp)]['dnabert_before']
            data_atac['dnabert_after'][i] = data_dnabert[(data_dnabert['POS']==cpg)&(data_dnabert['TSS']==snp)]['dnabert_after']
            data_atac['atac_between_ori'][i] = data_atac_ori[(data_atac_ori['POS']==cpg)&(data_atac_ori['TSS']==snp)]['atac_between']
        data_atac.to_pickle(merge_path + chr + '_' + model + '.dataset')

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

    train_dnabert = train_data[['RSID','RHO','A1','A2','CHR','TSS','POS','label','distance','seq_len','atac_between','dnabert_before','dnabert_after']]
    valid_dnabert = valid_data[['RSID','RHO','A1','A2','CHR','TSS','POS','label','distance','seq_len','atac_between','dnabert_before','dnabert_after']]
    test_dnabert = test_data[['RSID','RHO','A1','A2','CHR','TSS','POS','label','distance','seq_len','atac_between','dnabert_before','dnabert_after']]
    train_dnabert.to_pickle(output_path1 + model + '_train.dataset')
    valid_dnabert.to_pickle(output_path1 + model + '_valid.dataset')
    test_dnabert.to_pickle(output_path1 + model + '_test.dataset')