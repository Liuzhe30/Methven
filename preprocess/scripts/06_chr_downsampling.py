import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import json
import warnings
warnings.filterwarnings('ignore')

balanced_path = '../../datasets/balanced/'
downsample_path = '../../datasets/chr_downsampled/'

# scanning for min sample number 
count_dict = {'small':50_000,'middle':50_000,'large':50_000,'huge':50_000}
for i in range(22):
    chr = 'chr' + str(i+1)
    balanced_table = pd.read_pickle(balanced_path + chr + '_small.dataset')
    if(balanced_table.shape[0] < count_dict['small']):
        count_dict['small'] = balanced_table.shape[0]
    balanced_table = pd.read_pickle(balanced_path + chr + '_middle.dataset')
    if(balanced_table.shape[0] < count_dict['middle']):
        count_dict['middle'] = balanced_table.shape[0]
    balanced_table = pd.read_pickle(balanced_path + chr + '_large.dataset')
    if(balanced_table.shape[0] < count_dict['large']):
        count_dict['large'] = balanced_table.shape[0]
    balanced_table = pd.read_pickle(balanced_path + chr + '_huge.dataset')
    if(balanced_table.shape[0] < count_dict['huge']):
        count_dict['huge'] = balanced_table.shape[0]
print(count_dict) # {'small': 5024, 'middle': 6400, 'large': 9170, 'huge': 7896}

for i in range(22):
    chr = 'chr' + str(i+1)

    balanced_table = pd.read_pickle(balanced_path + chr + '_small.dataset')
    pos_table = balanced_table[balanced_table['label'] == 1].reset_index(drop=True)
    neg_table = balanced_table[balanced_table['label'] == 0].reset_index(drop=True)
    sample_num = int(count_dict['small']/6)
    pos_new_table = pos_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    neg_new_table = neg_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    final_table = pd.concat([pos_new_table,neg_new_table])
    final_table = final_table.sample(frac=1,replace=False).reset_index(drop=True)
    final_table.to_pickle(downsample_path + chr + '_small.dataset')

    balanced_table = pd.read_pickle(balanced_path + chr + '_middle.dataset')
    pos_table = balanced_table[balanced_table['label'] == 1].reset_index(drop=True)
    neg_table = balanced_table[balanced_table['label'] == 0].reset_index(drop=True)
    sample_num = int(count_dict['middle']/6)
    pos_new_table = pos_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    neg_new_table = neg_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    final_table = pd.concat([pos_new_table,neg_new_table])
    final_table = final_table.sample(frac=1,replace=False).reset_index(drop=True)
    final_table.to_pickle(downsample_path + chr + '_middle.dataset')

    balanced_table = pd.read_pickle(balanced_path + chr + '_large.dataset')
    pos_table = balanced_table[balanced_table['label'] == 1].reset_index(drop=True)
    neg_table = balanced_table[balanced_table['label'] == 0].reset_index(drop=True)
    sample_num = int(count_dict['large']/8)
    pos_new_table = pos_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    neg_new_table = neg_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    final_table = pd.concat([pos_new_table,neg_new_table])
    final_table = final_table.sample(frac=1,replace=False).reset_index(drop=True)
    final_table.to_pickle(downsample_path + chr + '_large.dataset')

    balanced_table = pd.read_pickle(balanced_path + chr + '_huge.dataset')
    pos_table = balanced_table[balanced_table['label'] == 1].reset_index(drop=True)
    neg_table = balanced_table[balanced_table['label'] == 0].reset_index(drop=True)
    sample_num = int(count_dict['huge']/14)
    pos_new_table = pos_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    neg_new_table = neg_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    final_table = pd.concat([pos_new_table,neg_new_table])
    final_table = final_table.sample(frac=1,replace=False).reset_index(drop=True)
    final_table.to_pickle(downsample_path + chr + '_huge.dataset')