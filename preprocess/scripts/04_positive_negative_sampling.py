import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import json
import warnings
warnings.filterwarnings('ignore')

cross_path = '../../datasets/middlefile/meQTL_annotation_CpG_SNP_chr/'
cpg_path = '../../datasets/middlefile/clean_epic/'
snp_dict_path = '../../datasets/middlefile/downsample_SNP.json'
save_path = '../../datasets/labeled/'

with open(snp_dict_path) as json_file:
    snp_dict = json.load(json_file)

'''
cross_table:
          CpG          SNP   Beta Ref Alt CHR   CpG_POS   SNP_POS  distance
0  cg11913416    rs1262461 -1.179   A   G   1  13582516  13581610       906
1  cg25722041  rs114812780  1.233   G   T   1   8563413   8468794     94619
2  cg25722041   rs12568293  1.233   T   C   1   8563413   8480852     82561
3  cg25722041   rs12117910  1.233   A   C   1   8563413   8496893     66520
4  cg25722041   rs12403339  1.233   G   A   1   8563413   8498232     65181

cpg_table:
       IlmnID CHR_hg38  Start_hg38   End_hg38
0  cg23229610     chr1     6781064    6781066
1  cg25458538     chr1   120548382  120548384
2  cg04118974     chr1   109816912  109816914
3  cg07659892     chr1   234523600  234523602
4  cg11993619     chr1   111743952  111743954
'''

for i in tqdm(range(22)):
    chr = 'chr' + str(i+1)
    cross_table = pd.read_pickle(cross_path + chr + '.pkl')
    cpg_table = pd.read_pickle(cpg_path + 'cpg_' + chr + '.pkl')
    cpg_table[['Start_hg38']] = cpg_table[['Start_hg38']].astype(int)

    # small model
    snp_list = snp_dict[chr]['small']
    dataset = pd.DataFrame(columns=['CpG', 'SNP', 'Beta', 'Ref', 'Alt', 'CHR', 'CpG_POS', 'SNP_POS','distance', 'label'])
    for snp in tqdm(snp_list):
        len_sample = cross_table[cross_table['distance'] <= 1_000].reset_index(drop=True)
        pos_sample = len_sample[len_sample['SNP'] == snp].reset_index(drop=True)
        snp_position = int(pos_sample['SNP_POS'][0])
        pos_cpg_list = pos_sample['CpG'].tolist()
        max_range = 1_000
        cpg_range = {'min':snp_position - max_range,'max':snp_position + max_range}
        candidate = cpg_table[(cpg_table['Start_hg38'] > cpg_range['min']) & (cpg_table['Start_hg38'] < cpg_range['max'])]
        for pos in pos_cpg_list:
            candidate_ = candidate[~(candidate['IlmnID'] == pos)].reset_index(drop=True)
        for j in range(pos_sample.shape[0]):
            dataset = dataset._append([{'CpG':pos_sample['CpG'][j], 'SNP':pos_sample['SNP'][0], 'Beta':pos_sample['Beta'][j], 'Ref':pos_sample['Ref'][0], 
                                    'Alt':pos_sample['Alt'][0], 'CHR':pos_sample['CHR'][0], 'CpG_POS':int(pos_sample['CpG_POS'][j]), 
                                    'SNP_POS':int(pos_sample['SNP_POS'][0]), 'distance':pos_sample['distance'][j],
                                    'label':1}], ignore_index=True)
        for j in range(candidate_.shape[0]):
            new_cpg_pos = int(candidate_['Start_hg38'][j])
            new_snp_pos = int(pos_sample['SNP_POS'][0])
            new_distance = np.abs(new_cpg_pos - new_snp_pos)
            dataset = dataset._append([{'CpG':candidate_['IlmnID'][j], 'SNP':pos_sample['SNP'][0], 'Beta':0, 'Ref':pos_sample['Ref'][0], 
                                    'Alt':pos_sample['Alt'][0], 'CHR':pos_sample['CHR'][0], 'CpG_POS':new_cpg_pos, 
                                    'SNP_POS':new_snp_pos, 'distance':new_distance,
                                    'label':0}], ignore_index=True)
    
    dataset.to_pickle(save_path + chr + '_small.dataset')
    print(dataset['label'].value_counts())
    #print(dataset.head())
    '''
                CpG         SNP    Beta Ref Alt CHR    CpG_POS    SNP_POS distance label
    0  cg19637330   rs4920348 -0.7066   C   T   1   18784427   18785300      873     1
    1  cg26220594   rs4920348 -0.6985   C   T   1   18784483   18785300      817     1
    2  cg14856563   rs4920348  0.0000   C   T   1   18784594   18785300      706     0
    3  cg19637330   rs4920348  0.0000   C   T   1   18784427   18785300      873     0
    4  cg21085632  rs12121609  0.3080   C   T   1  160101291  160101153      138     1
    '''

    # middle model
    snp_list = snp_dict[chr]['middle']
    dataset = pd.DataFrame(columns=['CpG', 'SNP', 'Beta', 'Ref', 'Alt', 'CHR', 'CpG_POS', 'SNP_POS','distance', 'label'])
    for snp in tqdm(snp_list):
        len_sample = cross_table[(cross_table['distance'] > 1_000) & (cross_table['distance'] <= 10_000)].reset_index(drop=True)
        pos_sample = len_sample[len_sample['SNP'] == snp].reset_index(drop=True)
        snp_position = int(pos_sample['SNP_POS'][0])
        pos_cpg_list = pos_sample['CpG'].tolist()
        max_range = 10_000
        min_range = 1_001
        cpg_range = {'min1':snp_position - max_range, 'max1':snp_position - min_range,
                     'min2':snp_position + min_range, 'max2':snp_position + max_range}
        candidate = cpg_table[(cpg_table['Start_hg38'] > cpg_range['min1']) & (cpg_table['Start_hg38'] < cpg_range['max1']) |
                               (cpg_table['Start_hg38'] > cpg_range['min2']) & (cpg_table['Start_hg38'] < cpg_range['max2'])]
        for pos in pos_cpg_list:
            candidate_ = candidate[~(candidate['IlmnID'] == pos)].reset_index(drop=True)
        for j in range(pos_sample.shape[0]):
            dataset = dataset._append([{'CpG':pos_sample['CpG'][j], 'SNP':pos_sample['SNP'][0], 'Beta':pos_sample['Beta'][j], 'Ref':pos_sample['Ref'][0], 
                                    'Alt':pos_sample['Alt'][0], 'CHR':pos_sample['CHR'][0], 'CpG_POS':int(pos_sample['CpG_POS'][j]), 
                                    'SNP_POS':int(pos_sample['SNP_POS'][0]), 'distance':pos_sample['distance'][j],
                                    'label':1}], ignore_index=True)
        for j in range(candidate_.shape[0]):
            new_cpg_pos = int(candidate_['Start_hg38'][j])
            new_snp_pos = int(pos_sample['SNP_POS'][0])
            new_distance = np.abs(new_cpg_pos - new_snp_pos)
            dataset = dataset._append([{'CpG':candidate_['IlmnID'][j], 'SNP':pos_sample['SNP'][0], 'Beta':0, 'Ref':pos_sample['Ref'][0], 
                                    'Alt':pos_sample['Alt'][0], 'CHR':pos_sample['CHR'][0], 'CpG_POS':new_cpg_pos, 
                                    'SNP_POS':new_snp_pos, 'distance':new_distance,
                                    'label':0}], ignore_index=True)

    dataset.to_pickle(save_path + chr + '_middle.dataset')
    print(dataset['label'].value_counts())
    #print(dataset.head())
    '''
             CpG         SNP    Beta Ref Alt CHR   CpG_POS   SNP_POS distance label                                                                                                                                                                                                                               | 20/2000 [00:01<01:41, 19.60it/s] 
    0  cg18376306  rs12044025  0.4642   A   C   1  75788822  75784030     4792     1
    1  cg00170796  rs12044025 -0.6113   A   C   1  75785851  75784030     1821     1
    2  cg05598546  rs12044025 -0.4703   A   C   1  75786080  75784030     2050     1
    3  cg13421759  rs12044025 -0.4591   A   C   1  75786066  75784030     2036     1
    4  cg07774299  rs12044025  0.3655   A   C   1  75787759  75784030     3729     1

              CpG         SNP Beta Ref Alt CHR   CpG_POS   SNP_POS distance label
    0  cg13250850  rs12044025    0   A   C   1  75785048  75784030     1018     0
    1  cg16429999  rs12044025    0   A   C   1  75786129  75784030     2099     0
    2  cg05598546  rs12044025    0   A   C   1  75786080  75784030     2050     0
    3  cg11245233  rs12044025    0   A   C   1  75786213  75784030     2183     0
    4  cg04940515  rs12044025    0   A   C   1  75786179  75784030     2149     0
    '''

    # large model
    snp_list = snp_dict[chr]['large']
    dataset = pd.DataFrame(columns=['CpG', 'SNP', 'Beta', 'Ref', 'Alt', 'CHR', 'CpG_POS', 'SNP_POS','distance', 'label'])
    for snp in tqdm(snp_list):
        len_sample = cross_table[(cross_table['distance'] > 10_000) & (cross_table['distance'] <= 100_000)].reset_index(drop=True)
        pos_sample = len_sample[len_sample['SNP'] == snp].reset_index(drop=True)
        snp_position = int(pos_sample['SNP_POS'][0])
        pos_cpg_list = pos_sample['CpG'].tolist()
        max_range = 100_000
        min_range = 10_001
        cpg_range = {'min1':snp_position - max_range, 'max1':snp_position - min_range,
                     'min2':snp_position + min_range, 'max2':snp_position + max_range}
        candidate = cpg_table[(cpg_table['Start_hg38'] > cpg_range['min1']) & (cpg_table['Start_hg38'] < cpg_range['max1']) |
                               (cpg_table['Start_hg38'] > cpg_range['min2']) & (cpg_table['Start_hg38'] < cpg_range['max2'])]
        for pos in pos_cpg_list:
            candidate_ = candidate[~(candidate['IlmnID'] == pos)].reset_index(drop=True)

        for j in range(pos_sample.shape[0]):
            dataset = dataset._append([{'CpG':pos_sample['CpG'][j], 'SNP':pos_sample['SNP'][0], 'Beta':pos_sample['Beta'][j], 'Ref':pos_sample['Ref'][0], 
                                    'Alt':pos_sample['Alt'][0], 'CHR':pos_sample['CHR'][0], 'CpG_POS':int(pos_sample['CpG_POS'][j]), 
                                    'SNP_POS':int(pos_sample['SNP_POS'][0]), 'distance':pos_sample['distance'][j],
                                    'label':1}], ignore_index=True)
        for j in range(candidate_.shape[0]):
            new_cpg_pos = int(candidate_['Start_hg38'][j])
            new_snp_pos = int(pos_sample['SNP_POS'][0])
            new_distance = np.abs(new_cpg_pos - new_snp_pos)
            dataset = dataset._append([{'CpG':candidate_['IlmnID'][j], 'SNP':pos_sample['SNP'][0], 'Beta':0, 'Ref':pos_sample['Ref'][0], 
                                    'Alt':pos_sample['Alt'][0], 'CHR':pos_sample['CHR'][0], 'CpG_POS':new_cpg_pos, 
                                    'SNP_POS':new_snp_pos, 'distance':new_distance,
                                    'label':0}], ignore_index=True)

    dataset.to_pickle(save_path + chr + '_large.dataset')
    print(dataset['label'].value_counts())
    #print(dataset.head())
    '''
             CpG         SNP    Beta Ref Alt CHR   CpG_POS   SNP_POS distance label                                                                                                                                                                                                                                | 6/2000 [00:01<07:03,  4.71it/s] 
    0  cg10676837  rs10158705 -0.6795   G   A   1  19452577  19424086    28491     1
    1  cg17081867  rs10158705 -0.4874   G   A   1  19441601  19424086    17515     1
    2  cg08558153  rs10158705 -0.2824   G   A   1  19446335  19424086    22249     1
    3  cg25159064  rs10158705 -0.2571   G   A   1  19390988  19424086    33098     1
    4  cg01832549  rs10158705  0.2954   G   A   1  19448494  19424086    24408     1

             CpG         SNP Beta Ref Alt CHR   CpG_POS   SNP_POS distance label                                                                                                                                                                                                                                   | 6/2000 [00:01<07:13,  4.60it/s] 
    0  cg15100209  rs10158705    0   G   A   1  19343951  19424086    80135     0
    1  cg11888066  rs10158705    0   G   A   1  19484750  19424086    60664     0
    2  cg15624109  rs10158705    0   G   A   1  19489667  19424086    65581     0
    3  cg25159064  rs10158705    0   G   A   1  19390988  19424086    33098     0
    4  cg03076210  rs10158705    0   G   A   1  19485209  19424086    61123     0
    '''
    
    # huge model
    snp_list = snp_dict[chr]['huge']
    dataset = pd.DataFrame(columns=['CpG', 'SNP', 'Beta', 'Ref', 'Alt', 'CHR', 'CpG_POS', 'SNP_POS','distance', 'label'])
    for snp in tqdm(snp_list):
        len_sample = cross_table[(cross_table['distance'] > 100_000) & (cross_table['distance'] <= 1_000_000)].reset_index(drop=True)
        pos_sample = len_sample[len_sample['SNP'] == snp].reset_index(drop=True)
        snp_position = int(pos_sample['SNP_POS'][0])
        pos_cpg_list = pos_sample['CpG'].tolist()
        max_range = 1_000_000
        min_range = 100_001
        cpg_range = {'min1':snp_position - max_range, 'max1':snp_position - min_range,
                     'min2':snp_position + min_range, 'max2':snp_position + max_range}
        candidate = cpg_table[(cpg_table['Start_hg38'] > cpg_range['min1']) & (cpg_table['Start_hg38'] < cpg_range['max1']) |
                               (cpg_table['Start_hg38'] > cpg_range['min2']) & (cpg_table['Start_hg38'] < cpg_range['max2'])]
        for pos in pos_cpg_list:
            candidate_ = candidate[~(candidate['IlmnID'] == pos)].reset_index(drop=True)

        for j in range(pos_sample.shape[0]):
            dataset = dataset._append([{'CpG':pos_sample['CpG'][j], 'SNP':pos_sample['SNP'][0], 'Beta':pos_sample['Beta'][j], 'Ref':pos_sample['Ref'][0], 
                                    'Alt':pos_sample['Alt'][0], 'CHR':pos_sample['CHR'][0], 'CpG_POS':int(pos_sample['CpG_POS'][j]), 
                                    'SNP_POS':int(pos_sample['SNP_POS'][0]), 'distance':pos_sample['distance'][j],
                                    'label':1}], ignore_index=True)
        for j in range(candidate_.shape[0]):
            new_cpg_pos = int(candidate_['Start_hg38'][j])
            new_snp_pos = int(pos_sample['SNP_POS'][0])
            new_distance = np.abs(new_cpg_pos - new_snp_pos)
            dataset = dataset._append([{'CpG':candidate_['IlmnID'][j], 'SNP':pos_sample['SNP'][0], 'Beta':0, 'Ref':pos_sample['Ref'][0], 
                                    'Alt':pos_sample['Alt'][0], 'CHR':pos_sample['CHR'][0], 'CpG_POS':new_cpg_pos, 
                                    'SNP_POS':new_snp_pos, 'distance':new_distance,
                                    'label':0}], ignore_index=True)

    dataset.to_pickle(save_path + chr + '_huge.dataset')
    print(dataset['label'].value_counts())