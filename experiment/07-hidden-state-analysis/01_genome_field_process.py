import pandas as pd

file_path = 'genome_all_field.tsv'
output_path = 'genome_field/'
data = pd.read_csv(file_path,sep='\t')
print(data.head())
'''
   #bin chrom  chromStart   chromEnd               name  score strand  thickStart   thickEnd      itemRgb
0     0  chr1    67001612   67113812  13_Heterochrom/lo      0      .    67001612   67113812  245,245,245
1     0  chr1   201324577  201326777       12_Repressed      0      .   201324577  201326777  127,127,127
2     1  chr1     8381013    8408813  13_Heterochrom/lo      0      .     8381013    8408813  245,245,245
3     1  chr1    16777013   16777413  5_Strong_Enhancer      0      .    16777013   16777413    250,202,0
4     1  chr1    25162013   25166013        11_Weak_Txn      0      .    25162013   25166013  153,255,10
'''

# list all types
print(data['name'].unique())
'''
['13_Heterochrom/lo' '12_Repressed' '5_Strong_Enhancer' '11_Weak_Txn' '6_Weak_Enhancer' '7_Weak_Enhancer' '10_Txn_Elongation'
 '14_Repetitive/CNV' '9_Txn_Transition' '1_Active_Promoter' '2_Weak_Promoter' '4_Strong_Enhancer' '8_Insulator' '3_Poised_Promoter'
 '15_Repetitive/CNV']
'''

# select cpg associated regions
type_dict = {
    '1_Active_Promoter':'active_promoter',
    '4_Strong_Enhancer':'strong_enhancer',
    '5_Strong_Enhancer':'strong_enhancer',
    '9_Txn_Transition':'txn_transition',
    '10_Txn_Elongation':'txn_elongation',
    '8_Insulator':'insulator',
    '13_Heterochrom/lo':'heterochrom_lo',
    '12_Repressed':'repressed',
    '14_Repetitive/CNV':'repetitive_CNV',
    '15_Repetitive/CNV':'repetitive_CNV'
}


# filter annotations
chr_list = data['chrom'].unique()
print(chr_list)
for c in chr_list:
    chr_data = data[data['chrom'] == c]
    chr_data = chr_data.reset_index(drop=True)
    print(chr_data.head())
    new_chr_data = pd.DataFrame()
    for i in range(chr_data.shape[0]):
        type = chr_data['name'][i]
        if(type in type_dict.keys()):
            new_chr_data = new_chr_data._append({'chrom':chr_data['chrom'][i],'chromStart':chr_data['chromStart'][i],
                                                    'chromEnd':chr_data['chromEnd'][i],'type':type_dict[type]},ignore_index=True)
    print(new_chr_data.head())
    new_chr_data.to_pickle(output_path + c + '_map.pkl')

'''
chrom  chromStart   chromEnd             type
0  chr1    67001612   67113812   heterochrom_lo
1  chr1   201324577  201326777        repressed
2  chr1     8381013    8408813   heterochrom_lo
3  chr1    16777013   16777413  strong enhancer
4  chr1    49401813   50512613   heterochrom_lo
'''

# get_annotation(chr):
import json
type_list = ['active_promoter','strong_enhancer','txn_transition',
    'txn_elongation','insulator','heterochrom_lo','repressed','repetitive_CNV']
for i in range(22):
    chr = 'chr' + str(i+1)
    annotation_file = pd.read_pickle(output_path + chr + '_map.pkl')
    annotation_dict = {}
    '''
    chrom  chromStart   chromEnd             type
    0  chr1    67001612   67113812   heterochrom_lo
    1  chr1   201324577  201326777        repressed
    2  chr1     8381013    8408813   heterochrom_lo
    3  chr1    16777013   16777413  strong enhancer
    4  chr1    49401813   50512613   heterochrom_lo
    '''
    for type in type_list:
        annotation_dict[type] = []
    for type in type_list:
        filtered = annotation_file[annotation_file['type']==type].reset_index(drop=True)
        for j in range(filtered.shape[0]):
            annotation_dict[type].append([int(filtered['chromStart'][j]),int(filtered['chromEnd'][j])])
    
    with open(output_path + chr + '_range.json','w') as f:
        json.dump(annotation_dict,f,indent=4)
