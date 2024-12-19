import os
import pandas as pd
import numpy as np
import subprocess
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

meqtl_path = '../../datasets/raw/mQTL_significant_pairs_retina.txt'
cpg_annotation_path = '../../datasets/middlefile/clean_450k/'
save_path = '../../datasets/middlefile/retina_meQTL_annotation/'

main_table = pd.read_csv(meqtl_path, sep='\t')
main_table = main_table.sample(frac=1/20, random_state=42).reset_index(drop=True)
main_table[['snp', 'SNP_POS', 'Ref', 'Alt']] = main_table['Variant'].str.split(':', expand=True)
main_table = main_table[(main_table['Ref'].str.len() == 1) & (main_table['Alt'].str.len() == 1)]

print(main_table.head())

for i in range(1, 23):  
    chr_num = 'chr' + str(i)
    chr_table = main_table[(main_table['CpG_Chr'] == 'chr' + str(i)) & (main_table['SNP_Chr'] == 'chr' + str(i))]
    new_table = pd.DataFrame(columns=['CpG', 'SNP', 'Beta', 'Ref', 'Alt', 'CHR', 'CpG_POS', 'SNP_POS', 'label'])

    for idx in tqdm(range(chr_table.shape[0])):
        try:
            cpg = chr_table.iloc[idx]['Phenotype']
            snp = chr_table.iloc[idx]['snp'] 
            snp_pos = chr_table.iloc[idx]['SNP_POS']
            ref = chr_table.iloc[idx]['Ref']
            alt = chr_table.iloc[idx]['Alt']
            beta = chr_table.iloc[idx]['Forward_slope']
            chr_label = chr_table.iloc[idx]['SNP_Chr']
            cpg_pos = chr_table.iloc[idx]['CpG_Start']
            snp_pos = chr_table.iloc[idx]['SNP_Start']

            label = 1 if beta > 0 else 0

            formatted_snp = f"{chr_table.iloc[idx]['SNP_Chr']}:{snp_pos}_{ref}_{alt}"

            new_row = {
                'CpG': cpg,
                'SNP': formatted_snp,
                'Beta': beta,
                'Ref': ref,
                'Alt': alt,
                'CHR': chr_label,
                'CpG_POS': cpg_pos,
                'SNP_POS': snp_pos,
                'label': label
            }
            new_table = new_table._append(new_row, ignore_index=True)

        except Exception as e:
            print(f"Error processing index {idx}: {e}")
            continue

    new_table.to_pickle(save_path + chr_num + '.pkl')
    print(f"File saved: {save_path + chr_num}.pkl")
    print(new_table.head())