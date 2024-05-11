# mapping ATAC-seq
import pandas as pd
import pyBigWig
import numpy as np
import math
pd.set_option('display.max_columns', None)

atac_path = '../../datasets/middlefile/kidney_atac.bw'
output_path = '../../datasets/atac_mapping/'
file_path = '../../datasets/seq_mapping_post/'

model_size = {'small':1_000,'middle':10_000,'large':100_000,'huge':1_000_000}
bw = pyBigWig.open(atac_path)
for i in range(22):
    chr = 'chr' + str(i+1)
    for model in model_size.keys():
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        max_range = model_size[model]
        data['atac_between'] = 0
        data['atac_between'] = data['atac_between'].astype('object')
        for i in range(len(data)):
            cpg_pos = int(data['CpG_POS'].values[i])
            snp_pos = int(data['SNP_POS'].values[i])
            seq_len = int(data['seq_len'].values[i])
            
            atac_between = bw.values(chr, cpg_pos - max_range - 1, cpg_pos + max_range)
            atac_between = [0 if math.isnan(x) else x for x in atac_between]
            data.at[i, 'atac_between'] = atac_between
        data.to_pickle(output_path + chr + '_' + model + '.dataset')