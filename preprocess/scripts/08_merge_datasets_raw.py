# merge dnabert embedding
import pandas as pd
import numpy as np
from sklearn.utils import shuffle
pd.set_option('display.max_columns', None)
import warnings
warnings.filterwarnings('ignore')

atac_path = '../../datasets/atac_mapping_post/'
output_path = '../../datasets/final/raw/'

model_size = {'small':10_000,'large':100_000}

# merge datasets
for model in model_size.keys():
    model_merged = pd.DataFrame()
    for i in range(22):
        chr = 'chr' + str(i+1)
        data = pd.read_pickle(atac_path + chr + '_' + model + '.dataset')
        model_merged = pd.concat([model_merged, data])
    model_merged = shuffle(model_merged)
    model_merged = model_merged.reset_index(drop=True)
    train_data = model_merged[0:int(0.8*len(model_merged))].reset_index(drop=True)
    valid_data = model_merged[int(0.8*len(model_merged)):int(0.9*len(model_merged))].reset_index(drop=True)
    test_data = model_merged[int(0.9*len(model_merged)):].reset_index(drop=True)
    train_data.to_pickle(output_path + model + '_train.dataset')
    valid_data.to_pickle(output_path + model + '_valid.dataset')
    test_data.to_pickle(output_path + model + '_test.dataset')