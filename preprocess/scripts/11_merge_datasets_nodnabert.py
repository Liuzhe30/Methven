# generate final datasets: onehot
# training: chr1-9,chr14-22
# validation: chr12-13
# test: chr10-11
import pandas as pd
pd.set_option('display.max_columns', None)

file_path = '../../datasets/atac_cut/'
output_path = '../../datasets/final/onehot/'

model_size = {'small':1_000,'middle':10_000,'large':100_000,'huge':1_000_000}

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