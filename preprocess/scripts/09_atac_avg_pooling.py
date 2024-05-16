import pandas as pd
import numpy as np
pd.set_option('display.max_columns', None)

file_path = '../../datasets/atac_mapping_post/'
output_path = '../../datasets/atac_cut/'

model_cutting = {'small':1,'middle':19,'large':199,'huge':1_999}
model_size = {'small':1_000,'middle':10_000,'large':100_000,'huge':1_000_000}

for i in range(22):
    chr = 'chr' + str(i+1)
    for model in model_size.keys():
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        for i in range(len(data)):
            atac = np.array(data['atac_between'][i])

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

            data['atac_between'][i] = new_atac_list
        data.to_pickle(output_path + chr + '_' + model + '.dataset')


