# mapping DNA sequences, centered on CpG
from sklearn.utils import shuffle
import pandas as pd
pd.set_option('display.max_columns', None)

cell_list = ['B_Naive','CD4_memory','CD8_memory']
for cell in cell_list:

    fasta_path = '../../datasets/raw/reference_genome_hg38/'
    file_path = '../../datasets/sampled_eQTL/' + cell 
    output_path = '../../datasets/seq_mapping_eQTL/' + cell 

    model_size = {'small':10_000,'large':100_000}
    for i in range(20):
        chr = 'chr' + str(i+1)
        for model in model_size.keys():
            data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
            data['seq_before'] = ''
            data['seq_after'] = ''
            data['seq_len'] = 0
            max_range = model_size[model]
            for i in range(len(data)):
                after_mutation = data['A2'].values[i]
                cpg_pos = int(data['POS'].values[i])
                snp_pos = int(data['TSS'].values[i])

                with open(fasta_path + chr + '.fasta') as fa:
                    line = fa.readline()
                    seq_before = line[cpg_pos - max_range - 1:cpg_pos + max_range]
                    seq_len = len(seq_before)
                    seq_after = line[cpg_pos - max_range - 1:snp_pos - 1] + after_mutation + line[snp_pos:cpg_pos + max_range] # replace A1

                data.loc[i, 'seq_before'] = seq_before
                data.loc[i, 'seq_after'] = seq_after
                data.loc[i, 'seq_len'] = seq_len
                print(data.head())
            data.to_pickle(output_path + chr + '_' + model + '.dataset')

    # post-process
    file_path = '../../datasets/seq_mapping_eQTL/' + cell 
    output_path = '../../datasets/seq_mapping_post_eQTL/' + cell 

    for i in range(20):
        chr = 'chr' + str(i+1)
        for model in model_size.keys():
            data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
            data_check = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
            max_range = model_size[model]
            for i in range(len(data_check)):
                seq_len = data_check['seq_len'][i]
                cpg = data_check['TSS'][i]
                snp = data_check['RSID'][i]
                if(seq_len != max_range * 2 + 1):
                    data = data.drop(data[(data['TSS']==cpg)&(data['RSID']==snp)].index)
            data = data.reset_index(drop=True)
            data.to_pickle(output_path + chr + '_' + model + '.dataset')