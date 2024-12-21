# mapping ATAC-seq
import pandas as pd
import pyBigWig
import numpy as np
import math
pd.set_option('display.max_columns', None)

cell_list = ['B_Naive','CD4_memory','CD8_memory']
for cell in cell_list:

    atac_path = '../../datasets/raw/' + cell + '.bigWig'
    output_path = '../../datasets/atac_mapping_eQTL/' + cell 
    file_path = '../../datasets/seq_mapping_post_eQTL/' + cell 

    model_size = {'small':10_000,'large':100_000}
    bw = pyBigWig.open(atac_path)
    for i in range(22):
        chr = 'chr' + str(i+1)
        for model in model_size.keys():
            data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
            max_range = model_size[model]
            data['atac_between'] = 0
            data['atac_between'] = data['atac_between'].astype('object')
            for i in range(len(data)):
                cpg_pos = int(data['POS'].values[i])
                snp_pos = int(data['TSS'].values[i])
                seq_len = int(data['seq_len'].values[i])
                
                atac_between = bw.values(chr, cpg_pos - max_range - 1, cpg_pos + max_range)
                atac_between = [0 if math.isnan(x) else x for x in atac_between]
                data.at[i, 'atac_between'] = atac_between
            data.to_pickle(output_path + chr + '_' + model + '.dataset')

    # post-process

    def check(check_list):
        for elem in check_list:
            if(float(elem) == np.nan):
                return 0
            if(math.isnan(elem) == True):
                return 0
        return 1

    file_path = '../../datasets/atac_mapping_eQTL/' + cell 
    output_path = '../../datasets/atac_mapping_post_eQTL/' + cell 

    for i in range(22):
        chr = 'chr' + str(i+1)
        for model in model_size.keys():
            data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
            data_check = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
            for i in range(len(data_check)):
                atac_between = data_check['atac_between'][i]
                cpg = data_check['TSS'][i]
                snp = data_check['RSID'][i]
                if(check(atac_between)==0):
                    data = data.drop(data[(data['TSS']==cpg)&(data['RSID']==snp)].index)
            data = data.reset_index(drop=True)
            data.to_pickle(output_path + chr + '_' + model + '.dataset')

    # average pooling
    file_path = '../../datasets/atac_mapping_post_eQTL/' + cell 
    output_path = '../../datasets/atac_cut_eQTL/' + cell 

    model_cutting = {'small':19,'large':199}

    for c in range(22):
        chr = 'chr' + str(c+1)
        for model in model_size.keys():
            data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
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
            data.to_pickle(output_path + chr + '_' + model + '.dataset')