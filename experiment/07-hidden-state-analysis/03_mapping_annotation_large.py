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

# model 2: large model
n = 199 # position-wise cutting
half_range = 100_000
data = pd.read_pickle(dataset_path + 'large.dataset')
'''
          CpG               SNP      Beta    CHR    CpG_POS    SNP_POS                                       bigru_output  distance
0  cg06461408   17:73599603_C_T  6.548007  chr17   73561668   73599603  [[-0.1599046, -0.08208149, -0.10797326, -0.021...     37935
1  cg12797254   14:61764115_T_C  2.610837  chr14   61825961   61764115  [[-0.006501222, 0.031860273, -0.077744216, -0....     61846
2  cg03840920   16:12237834_G_A -3.787189  chr16   12321227   12237834  [[-0.20910944, 0.048690807, -0.13291685, -0.18...     83393
3  cg07569984   14:24843620_T_C -3.014296  chr14   24759069   24843620  [[-0.12150787, -0.012298661, -0.22109702, -0.0...     84551
4  cg13272119  10:131534075_A_G -5.001347  chr10  131460392  131534075  [[0.048231304, -0.008969577, -0.14642906, -0.2...     73683
'''
print(data['bigru_output'][0].shape) # (401, 128)

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

    np.save(output_path + type + '_large.npy',type_array)
    np.save(output_path + type + '_not_large.npy',not_type_array)