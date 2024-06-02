# generate 10-fold datasets, for the sake of partition consistency, data sets are no longer shuffled
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', None)
import warnings
warnings.filterwarnings('ignore')

file_path1 = '../../datasets/final/dnabert/'
file_path2 = '../../datasets/final/origin/'
output_path1 = '../../datasets/ten_fold/dnabert/'
output_path2 = '../../datasets/ten_fold/origin/'
clean_path = '../../data/ten_fold/'

model_size = {'small':10_000,'large':100_000}

for model in model_size.keys():
    train_data = pd.read_pickle(file_path1 + model + '_train.dataset')
    spl_df = np.array_split(train_data,10)
    for i in range(10):
        sub_df = spl_df[i].reset_index(drop=True)
        print(sub_df.shape)
        sub_df.to_pickle(output_path1 + model + '_split' + str(i+1) + '.dataset')
        clean_df = sub_df[['CpG','SNP','Beta','Ref','Alt','CHR','CpG_POS','SNP_POS']]
        clean_df.to_pickle(clean_path + model + '_split' + str(i+1) + '.dataset')
    train_data = pd.read_pickle(file_path2 + model + '_train.dataset')
    spl_df = np.array_split(train_data,10)
    for i in range(10):
        sub_df = spl_df[i].reset_index(drop=True)
        print(sub_df.shape)
        sub_df.to_pickle(output_path2 + model + '_split' + str(i+1) + '.dataset')