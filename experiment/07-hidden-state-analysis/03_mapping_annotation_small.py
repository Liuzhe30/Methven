import pandas as pd
import numpy as np
import json
from tqdm import tqdm

dataset_path = 'hidden_state/'
annotation_path = 'genome_field/'
output_path = 'results/'
type_list = ['active_promoter','strong_enhancer','txn_transition',
    'txn_elongation','insulator','heterochrom_lo','repressed','repetitive_CNV']

def generate_coordinate(n, half_range, cpg_pos,bigru_output):
    coord_dict = {}
    start_pos = cpg_pos - half_range - 1
    # first 250bp
    for i in range(start_pos, start_pos+250):
        coord_dict[i] = bigru_output[0].sum()
    # first n*500bp
    for j in range(n):
        for i in range(start_pos+250+500*j,start_pos+250+500*(j+1)):
            coord_dict[i] = bigru_output[j+1].sum()
    # mutation part
    for i in range(start_pos+250+500*n,start_pos+250+500*n+501):
        coord_dict[i] = bigru_output[n+1].sum()
    # second n*500bp
    for j in range(n):
        for i in range(start_pos+250+500*n+501+500*j,start_pos+250+500*n+501+500*(j+1)):
            coord_dict[i] = bigru_output[n+j+1].sum()
    # last n*500bp
    for i in range(start_pos+250+500*n+501+500*n, 250+500*n+501+500*n+250):
        coord_dict[i] = bigru_output[n+n+1+1].sum()

    return coord_dict

# model 1: small model
n = 19 # position-wise cutting
half_range = 10_000
data = pd.read_pickle(dataset_path + 'small.dataset')
'''
          CpG              SNP      Beta    CHR    CpG_POS    SNP_POS                                       bigru_output  distance
0  cg18621232   1:31675417_A_G  3.461477   chr1   31681696   31675417  [[0.12180035, -0.045385737, 0.12046003, 0.0148...      6279
1  cg22603971  15:38844106_G_A  3.578342  chr15   38850128   38844106  [[0.14448221, -0.032808937, 0.12122999, 0.0015...      6022
2  cg09447675  12:46877429_A_G  1.575897  chr12   46877466   46877429  [[0.021617461, 0.010252297, 0.15253885, -0.032...        37
3  cg05219430  8:122000893_G_A -2.366001   chr8  122005631  122000893  [[0.07814622, -0.08472212, 0.17048402, 0.00085...      4738
4  cg13033417   8:17694702_A_G  3.640821   chr8   17689804   17694702  [[-0.10019994, 0.2363034, 0.13702674, 0.174779...      4898
'''
print(data['bigru_output'][0].shape) # (41, 128)

count_dict = {}
for type in type_list:
    count_dict[type] = []
    count_dict['not_' + type] = []

for i in tqdm(range(data.shape[0])):
    chr = data['CHR'][i]
    cpg_pos = data['CpG_POS'][i]
    bigru_output = data['bigru_output'][i]
    coord_dict = generate_coordinate(n,half_range,cpg_pos,bigru_output)

    with open(annotation_path + chr + '_range.json','r',encoding='utf-8') as load_f:
        load_dict = json.load(load_f)
    for type in type_list:
        for position in coord_dict.keys(): # DNA positions
            flag = 0
            for item in load_dict[type]: # many ranges
                if(position >= item[0] and position <= item[1]):
                    flag = 1
                    count_dict[type].append(coord_dict[position])
                    pass
            if(flag == 0):
                count_dict['not_' + type].append(coord_dict[position])

for type in type_list:
    type_array = np.array(count_dict[type])
    not_type_array = np.array(count_dict['not_' + type])

    print(type)
    print(type_array.shape)
    print(not_type_array.shape)

    np.save(output_path + type + '_small.npy',type_array)
    np.save(output_path + type + '_not_small.npy',not_type_array)