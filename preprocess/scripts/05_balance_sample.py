import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import json
import warnings
warnings.filterwarnings('ignore')

labeled_path = '../../datasets/labeled/'
balanced_path = '../../datasets/balanced/'

for i in range(22):
    chr = 'chr' + str(i+1)
    # small table
    labeled_table = pd.read_pickle(labeled_path + chr + '_small.dataset')
    pos_table = labeled_table[labeled_table['label'] == 1].reset_index(drop=True)
    sample_num = pos_table.shape[0]
    neg_table = labeled_table[labeled_table['label'] == 0].reset_index(drop=True)
    neg_balance_table = neg_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    final_table = pd.concat([pos_table,neg_balance_table])
    final_table = final_table.sample(frac=1,replace=False).reset_index(drop=True)
    final_table.to_pickle(balanced_path + chr + '_small.dataset')
    #print(final_table.head())
    '''
          CpG         SNP    Beta Ref Alt CHR   CpG_POS   SNP_POS distance label
0  cg25824828  rs12170452  0.0000   A   G  22  39623110  39623768      658     0
1  cg03407441     rs81027 -0.2642   C   T  22  20956438  20957243      805     1
2  cg20254179    rs140573  0.4649   A   G  22  45056127  45056428      301     1
3  cg21716444   rs7290274  0.0000   T   C  22  37519242  37519908      666     0
4  cg07897659   rs9616213  0.0000   T   G  22  49923252  49924203      951     0
    '''

    # middle table
    labeled_table = pd.read_pickle(labeled_path + chr + '_middle.dataset')
    pos_table = labeled_table[labeled_table['label'] == 1].reset_index(drop=True)
    sample_num = pos_table.shape[0]
    neg_table = labeled_table[labeled_table['label'] == 0].reset_index(drop=True)
    neg_balance_table = neg_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    final_table = pd.concat([pos_table,neg_balance_table])
    final_table = final_table.sample(frac=1,replace=False).reset_index(drop=True)
    final_table.to_pickle(balanced_path + chr + '_middle.dataset')

    # large table
    labeled_table = pd.read_pickle(labeled_path + chr + '_large.dataset')
    pos_table = labeled_table[labeled_table['label'] == 1].reset_index(drop=True)
    sample_num = pos_table.shape[0]
    neg_table = labeled_table[labeled_table['label'] == 0].reset_index(drop=True)
    neg_balance_table = neg_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    final_table = pd.concat([pos_table,neg_balance_table])
    final_table = final_table.sample(frac=1,replace=False).reset_index(drop=True)
    final_table.to_pickle(balanced_path + chr + '_large.dataset')

    # huge table
    labeled_table = pd.read_pickle(labeled_path + chr + '_huge.dataset')
    pos_table = labeled_table[labeled_table['label'] == 1].reset_index(drop=True)
    sample_num = pos_table.shape[0]
    neg_table = labeled_table[labeled_table['label'] == 0].reset_index(drop=True)
    neg_balance_table = neg_table.sample(n=sample_num,replace=False).reset_index(drop=True)
    final_table = pd.concat([pos_table,neg_balance_table])
    final_table = final_table.sample(frac=1,replace=False).reset_index(drop=True)
    final_table.to_pickle(balanced_path + chr + '_huge.dataset')
