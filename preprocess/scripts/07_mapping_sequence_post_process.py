# post-process after mapping DNA sequences
import pandas as pd
pd.set_option('display.max_columns', None)

file_path = '../../datasets/seq_mapping/'
output_path = '../../datasets/seq_mapping_post/'

model_size = {'small':1_000,'middle':10_000,'large':100_000,'huge':1_000_000}
for i in range(22):
    chr = 'chr' + str(i+1)
    for model in model_size.keys():
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        data_check = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        max_range = model_size[model]
        for i in range(len(data_check)):
            seq_len = data_check['seq_len'][i]
            cpg = data_check['CpG'][i]
            snp = data_check['SNP'][i]
            if(seq_len != max_range * 2 + 1):
                data = data.drop(data[(data['CpG']==cpg)&(data['SNP']==snp)].index)
        data = data.reset_index(drop=True)
        data.to_pickle(output_path + chr + '_' + model + '.dataset')