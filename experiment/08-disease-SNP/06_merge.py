# merge dnabert embedding
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
import warnings
warnings.filterwarnings('ignore')

atac_path = 'data/atac_cut_t0/'
dnabert_path = 'data/dnabert_embedding/'
merge_path = 'data/final_t0/'

snp_list = ['rs2476601','rs3806624', 'rs7731626', 'rs2234067','rs2233424','rs947474','rs3824660','rs968567','rs3218251']
model_size = {'small':10_000,'large':100_000}

for snp in snp_list:
    for model in model_size.keys():
        data_atac = pd.read_pickle(atac_path + snp + '_' + model + '.dataset')
        data_dnabert = pd.read_pickle(dnabert_path + snp + '_' + model + '.dataset')
        data_atac['dnabert_before'] = 0
        data_atac['dnabert_before'] = data_atac['dnabert_before'].astype('object')
        data_atac['dnabert_after'] = 0
        data_atac['dnabert_after'] = data_atac['dnabert_after'].astype('object')

        for i in range(len(data_atac)):
            cpg = data_atac['CpG'][i]
            snp = data_atac['SNP'][i]
            data_atac['dnabert_before'][i] = data_dnabert[(data_dnabert['CpG']==cpg)&(data_dnabert['SNP']==snp)]['dnabert_before']
            data_atac['dnabert_after'][i] = data_dnabert[(data_dnabert['CpG']==cpg)&(data_dnabert['SNP']==snp)]['dnabert_after']
        data_atac.to_pickle(merge_path + snp + '_' + model + '.dataset')


atac_path = 'data/atac_cut_t24/'
dnabert_path = 'data/dnabert_embedding/'
merge_path = 'data/final_t24/'

for snp in snp_list:
    for model in model_size.keys():
        data_atac = pd.read_pickle(atac_path + snp + '_' + model + '.dataset')
        data_dnabert = pd.read_pickle(dnabert_path + snp + '_' + model + '.dataset')
        data_atac['dnabert_before'] = 0
        data_atac['dnabert_before'] = data_atac['dnabert_before'].astype('object')
        data_atac['dnabert_after'] = 0
        data_atac['dnabert_after'] = data_atac['dnabert_after'].astype('object')

        for i in range(len(data_atac)):
            cpg = data_atac['CpG'][i]
            snp = data_atac['SNP'][i]
            data_atac['dnabert_before'][i] = data_dnabert[(data_dnabert['CpG']==cpg)&(data_dnabert['SNP']==snp)]['dnabert_before']
            data_atac['dnabert_after'][i] = data_dnabert[(data_dnabert['CpG']==cpg)&(data_dnabert['SNP']==snp)]['dnabert_after']
        data_atac.to_pickle(merge_path + snp + '_' + model + '.dataset')