# generate clean datasets & fetch only sequences
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', None)

file_path = '../../datasets/merge_dnabert/'
output_path = '../../datasets/clean_final_datasets/'

model_size = {'small':1_000,'middle':10_000,'large':100_000,'huge':1_000_000}
for i in range(22):
    chr = 'chr' + str(i+1)
    for model in model_size.keys():
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        new_data = data.drop(columns=['seq_before','seq_after','atac_between','dnabert_before','dnabert_after','seq_len'])
        new_data.to_pickle(output_path + chr + '_' + model + '.dataset')
