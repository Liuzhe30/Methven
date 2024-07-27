import pandas as pd
import numpy as np

state_path = '../../model/middle_output_biGRU/'
dataset_path = '../../data/original_split/'
save_path = 'hidden_state/'

model_list = ['small','large']

'''
          CpG              SNP      Beta Ref Alt    CHR    CpG_POS    SNP_POS
0  cg18621232   1:31675417_A_G  3.461477   A   G   chr1   31681696   31675417
1  cg22603971  15:38844106_G_A  3.578342   G   A  chr15   38850128   38844106
2  cg09447675  12:46877429_A_G  1.575897   A   G  chr12   46877466   46877429
3  cg05219430  8:122000893_G_A -2.366001   G   A   chr8  122005631  122000893
4  cg13033417   8:17694702_A_G  3.640821   A   G   chr8   17689804   17694702
'''

for model in model_list:
    train_data = pd.read_pickle(dataset_path + model + '_train.dataset')
    state_array = np.load(state_path + model + '_output.npy')
    print(state_array.shape) # (15899, 41, 128)

    print(train_data.head())
    train_data = train_data[['CpG','SNP','Beta','CHR','CpG_POS','SNP_POS']]
    train_data['bigru_output'] = 0
    train_data['bigru_output'] = train_data['bigru_output'].astype('object')
    train_data['distance'] = 0

    for i in range(train_data.shape[0]):
        train_data['bigru_output'][i] = state_array[i]
        cpg_pos = train_data['CpG_POS'][i]
        snp_pos = train_data['SNP_POS'][i]
        distance = np.abs(int(cpg_pos) - int(snp_pos))
        train_data['distance'][i] = distance
    
    print(train_data.head())
    train_data.to_pickle(save_path + model + '.dataset')

'''
          CpG              SNP      Beta    CHR    CpG_POS    SNP_POS                                       bigru_output  distance
0  cg18621232   1:31675417_A_G  3.461477   chr1   31681696   31675417  [[0.12180035, -0.045385737, 0.12046003, 0.0148...      6279
1  cg22603971  15:38844106_G_A  3.578342  chr15   38850128   38844106  [[0.14448221, -0.032808937, 0.12122999, 0.0015...      6022
2  cg09447675  12:46877429_A_G  1.575897  chr12   46877466   46877429  [[0.021617461, 0.010252297, 0.15253885, -0.032...        37
3  cg05219430  8:122000893_G_A -2.366001   chr8  122005631  122000893  [[0.07814622, -0.08472212, 0.17048402, 0.00085...      4738
4  cg13033417   8:17694702_A_G  3.640821   chr8   17689804   17694702  [[-0.10019994, 0.2363034, 0.13702674, 0.174779...      4898
'''