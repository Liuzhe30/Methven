import os
import pandas as pd
import numpy as np
import subprocess
from tqdm import tqdm
import random
import json
import warnings
warnings.filterwarnings('ignore')

cell_list = ['B_Naive','CD4_memory','CD8_memory']
cell = cell_list[2]

save_path = '../../datasets/middlefile/eQTL_annotation/' + cell 
original_data = pd.read_pickle('../../datasets/raw/sceQTL/' + cell + '_hg38.pkl')
# split files
for i in range(22):
	chr = 'chr' + str(i+1)
	chr_table = original_data[original_data['CHR'] ==i+1].reset_index(drop=True)
	chr_table.to_pickle(save_path + chr + '.pkl')


cross_path = '../../datasets/middlefile/eQTL_annotation/' + cell 
save_path = '../../datasets/labeled_eQTL/' + cell 

'''
             RSID CHR    GENE A1 A2       RHO        TSS        POS label
0       rs9431818   1  GALNT2  C  T  0.107049  230067237  230282253     1
1       rs6541320   1    URB2  C  T  0.112008  229626246  230547657     1
2       rs4304537   1   TEX35  C  T -0.323398  178513108  178507419     0
3       rs4076972   1   TEX35  G  A -0.297878  178513108  178508297     0
4      rs10798610   1   TEX35  T  G -0.323339  178513108  178513849     0
'''

# add distance
for i in range(22):
	chr = 'chr' + str(i+1)
	cross_table = pd.read_pickle(cross_path + chr + '.pkl')
	cross_table['distance'] = 0
	for i in tqdm(range(cross_table.shape[0])):
		cpg_pos = cross_table['TSS'][i]
		snp_pos = cross_table['POS'][i]
		distance = np.abs(int(cpg_pos) - int(snp_pos))
		cross_table['distance'][i] = distance
		cross_table.to_pickle(cross_path + chr + '_distance.pkl')
	print(cross_table.head())   


# split by distance
cross_path = '../../datasets/middlefile/eQTL_annotation/' + cell 
save_path = '../../datasets/labeled_eQTL/' + cell 
for i in range(21):
	chr = 'chr' + str(i+1)
	cross_table = pd.read_pickle(cross_path + chr + '_distance.pkl')
    # small table
	len_sample = cross_table[cross_table['distance'] <= 10_000].reset_index(drop=True)
	len_sample.to_pickle(save_path  + chr + '_small.dataset')
	# large table
	len_sample = cross_table[(cross_table['distance'] > 10_000) & (cross_table['distance'] <= 100_000)].reset_index(drop=True)
	len_sample.to_pickle(save_path  + chr + '_large.dataset')

# down-sampling
cross_path = '../../datasets/labeled_eQTL/' + cell 
save_path = '../../datasets/sampled_eQTL/' + cell 
for i in range(21):
	chr = 'chr' + str(i+1)
	cross_table1 = pd.read_pickle(cross_path  + chr + '_small.dataset')
	cross_table2 = pd.read_pickle(cross_path  + chr + '_large.dataset')
	if(chr == 'chr6'):
		# small table
		min_sample = min(cross_table1['label'].value_counts()[0], cross_table1['label'].value_counts()[1])
		min_sample = int(min_sample/2)
		up_table = cross_table1[cross_table1['label']==1].sample(frac=1,replace=False).reset_index(drop=True)[0:min_sample]
		down_table = cross_table1[cross_table1['label']==0].sample(frac=1,replace=False).reset_index(drop=True)[0:min_sample]
		new_table = pd.concat([up_table,down_table]).reset_index(drop=True)
		#print(new_table['label'].value_counts())
		new_table.to_pickle(save_path  + chr + '_small.dataset')
		# large table
		min_sample = min(cross_table2['label'].value_counts()[0], cross_table2['label'].value_counts()[1])
		min_sample = int(min_sample/5)
		up_table = cross_table2[cross_table2['label']==1].sample(frac=1,replace=False).reset_index(drop=True)[0:min_sample]
		down_table = cross_table2[cross_table2['label']==0].sample(frac=1,replace=False).reset_index(drop=True)[0:min_sample]
		new_table = pd.concat([up_table,down_table]).reset_index(drop=True)
		#print(new_table['label'].value_counts())
		new_table.to_pickle(save_path  + chr + '_large.dataset')
	else:
		cross_table1.to_pickle(save_path  + chr + '_small.dataset')
		cross_table2.to_pickle(save_path  + chr + '_large.dataset')
