# mapping ATAC-seq
import pandas as pd
import pyBigWig
import numpy as np
import math
pd.set_option('display.max_columns', None)

atac_path = 'hg19_ATAC/hg19_t24h_atac.bw'
output_path = 'data/atac_mapping_t24/'
file_path = 'data/seq_mapping_post/'

snp_list = ['rs2476601','rs3806624', 'rs7731626', 'rs2234067','rs2233424','rs947474','rs3824660','rs968567','rs3218251']

model_size = {'small':10_000,'large':100_000}
bw = pyBigWig.open(atac_path)
for snp in snp_list:
    for model in model_size.keys():
        data = pd.read_pickle(file_path + snp + '_' + model + '.dataset')
        max_range = model_size[model]
        data['atac_between'] = 0
        data['atac_between'] = data['atac_between'].astype('object')
        for i in range(len(data)):
            chr = 'chr' + str(data['CHR'].values[i])
            cpg_pos = int(data['CpG_POS'].values[i])
            snp_pos = int(data['SNP_POS'].values[i])
            seq_len = int(data['seq_len'].values[i])
            
            atac_between = bw.values(chr, cpg_pos - max_range - 1, cpg_pos + max_range)
            atac_between = [0 if math.isnan(x) else x for x in atac_between]
            data.at[i, 'atac_between'] = atac_between
        data.to_pickle(output_path + snp + '_' + model + '.dataset')

# post-process

def check(check_list):
    for elem in check_list:
        if(float(elem) == np.nan):
            return 0
        if(math.isnan(elem) == True):
            return 0
    return 1

file_path = 'data/atac_mapping_t24/'
output_path = 'data/atac_mapping_post_t24/'

for snp in snp_list:
    print(snp)
    for model in model_size.keys():
        print(model)
        data = pd.read_pickle(file_path + snp + '_' + model + '.dataset')
        data_check = pd.read_pickle(file_path + snp + '_' + model + '.dataset')
        for i in range(len(data_check)):
            atac_between = data_check['atac_between'][i]
            cpg = data_check['CpG'][i]
            snp = data_check['SNP'][i]
            if(check(atac_between)==0):
                data = data.drop(data[(data['CpG']==cpg)&(data['SNP']==snp)].index)
        data = data.reset_index(drop=True)
        data.to_pickle(output_path + snp + '_' + model + '.dataset')
        print(data.shape)

# average pooling
file_path = 'data/atac_mapping_post_t24/'
output_path = 'data/atac_cut_t24/'

model_cutting = {'small':19,'large':199}

for snp in snp_list:
    for model in model_size.keys():
        data = pd.read_pickle(file_path + snp + '_' + model + '.dataset')
        for idx in range(len(data)):
            atac = np.array(data['atac_between'][idx])

            new_atac_list = []
            n = model_cutting[model]
            # first 250bp
            new_atac_list.append(np.mean(atac[0:250]))
            # first n*500bp
            for i in range(n):
                new_atac_list.append(np.mean(atac[250+i*500:250+(i+1)*500]))
            # mutation part
            new_atac_list.append(np.mean(atac[250+500*n:250+500*n+501]))
            # second n*500bp
            for i in range(n):
                new_atac_list.append(np.mean(atac[250+500*n+501+i*500:250+500*n+501+(i+1)*500]))
            # second 250bp
            new_atac_list.append(np.mean(atac[250+500*n+501+500*n:250+500*n+501+500*n+250]))

            data['atac_between'][idx] = new_atac_list
        data.to_pickle(output_path + snp + '_' + model + '.dataset')