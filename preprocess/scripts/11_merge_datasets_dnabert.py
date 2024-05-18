# merge dnabert embedding
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', None)
import warnings
warnings.filterwarnings('ignore')

file_path = '../../datasets/atac_cut/'
dnabert_path = '/lustre/home/acct-bmelgn/bmelgn-4/methven/datasets/dnabert_embedding/'
merge_path = '../../datasets/merge_dnabert/'
output_path = '../../datasets/final/dnabert/'

model_size = {'small':1_000,'middle':10_000,'large':100_000,'huge':1_000_000}

for i in range(22):
    chr = 'chr' + str(i+1)
    for model in model_size.keys():
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        data_check = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        data['dnabert_before'] = 0
        data['dnabert_before'] = data['dnabert_before'].astype('object')
        data['dnabert_after'] = 0
        data['dnabert_after'] = data['dnabert_after'].astype('object')
        for i in range(len(data_check)):
            cpg = data_check['CpG'][i]
            snp = data_check['SNP'][i]
            try:
                data['dnabert_before'][i] = np.load(dnabert_path + model + '/' + cpg + '_' + snp + '_' + model + '_' + chr + '_before.npy')
                data['dnabert_after'][i] = np.load(dnabert_path + model + '/' + cpg + '_' + snp + '_' + model + '_' + chr + '_after.npy')
            except FileNotFoundError:
                data = data.drop(data[(data['CpG']==cpg)&(data['SNP']==snp)].index)
                print(cpg,snp)
        data.to_pickle(merge_path + chr + '_' + model + '.dataset')

# merge
for model in model_size.keys():
    train_data = pd.DataFrame()
    valid_data = pd.DataFrame()
    test_data = pd.DataFrame()
    for i in range(1,10):
        chr = 'chr' + str(i)
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        train_data = pd.concat([train_data,data])
    for i in range(10,12):
        chr = 'chr' + str(i)
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        test_data = pd.concat([test_data,data])
    for i in range(12,14):
        chr = 'chr' + str(i)
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        valid_data = pd.concat([valid_data,data])
    for i in range(14,23):
        chr = 'chr' + str(i)
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        train_data = pd.concat([train_data,data])
    train_data = train_data.reset_index(drop=True)
    valid_data = valid_data.reset_index(drop=True)
    test_data = test_data.reset_index(drop=True)
    train_data.to_pickle(output_path + 'train_' + model + '.dataset')
    valid_data.to_pickle(output_path + 'valid_' + model + '.dataset')
    test_data.to_pickle(output_path + 'test_' + model + '.dataset')