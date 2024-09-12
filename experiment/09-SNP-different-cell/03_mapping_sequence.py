# mapping DNA sequences, centered on CpG
import pandas as pd
pd.set_option('display.max_columns', None)

fasta_path = '../../datasets/raw/reference_genome_hg19/'
file_path = 'data/snp_cpg/'
output_path = 'data/seq_mapping/'

snp_list = ['rs3806624', 'rs7731626', 'rs2234067','rs2233424','rs947474','rs3824660','rs968567','rs3218251']

'''
          CpG        SNP  Beta Ref Alt CHR    CpG_POS    SNP_POS  label  distance
0  cg12900790  rs2476601     0   A   G   1  114375675  114377568      0      1893
1  cg01439475  rs2476601     0   A   G   1  114372756  114377568      0      4812
'''

model_size = {'small':10_000,'large':100_000}
for snp in snp_list:
    for model in model_size.keys():
        data = pd.read_pickle(file_path + snp + '_' + model + '.pkl')
        data['seq_before'] = ''
        data['seq_after'] = ''
        data['seq_len'] = 0
        max_range = model_size[model]
        for i in range(len(data)):
            after_mutation = data['Alt'].values[i]
            cpg_pos = int(data['CpG_POS'].values[i])
            snp_pos = int(data['SNP_POS'].values[i])
            chr = 'chr' + str(data['CHR'].values[i])

            with open(fasta_path + chr + '.fasta') as fa:
                line = fa.readline()
                seq_before = line[cpg_pos - max_range - 1:cpg_pos + max_range]
                seq_len = len(seq_before)
                seq_after = line[cpg_pos - max_range - 1:snp_pos - 1] + after_mutation + line[snp_pos:cpg_pos + max_range] # replace A1

            data.loc[i, 'seq_before'] = seq_before
            data.loc[i, 'seq_after'] = seq_after
            data.loc[i, 'seq_len'] = seq_len
        print(data.head())
        data.to_pickle(output_path + snp + '_' + model + '.dataset')

'''
          CpG        SNP  Beta Ref Alt CHR    CpG_POS    SNP_POS  label  \
0  cg12900790  rs2476601     0   A   G   1  114375675  114377568      0
1  cg01439475  rs2476601     0   A   G   1  114372756  114377568      0

   distance                                         seq_before  \
0      1893  tatcatgtgttattaagtatatatttgtgaatatggtcaaccaaat...
1      4812  TTTTTTTCTAATTATAAGTGATTCTAATTTAGGCCCATAAccgggc...

                                           seq_after  seq_len
0  tatcatgtgttattaagtatatatttgtgaatatggtcaaccaaat...    20001
1  TTTTTTTCTAATTATAAGTGATTCTAATTTAGGCCCATAAccgggc...    20001
'''
'''
          CpG        SNP  Beta Ref Alt CHR    CpG_POS    SNP_POS  label  \
0  cg08266106  rs2476601     0   A   G   1  114430903  114377568      0
1  cg05825720  rs2476601     0   A   G   1  114301557  114377568      0
2  cg02066222  rs2476601     0   A   G   1  114354688  114377568      0
3  cg13572289  rs2476601     0   A   G   1  114447746  114377568      0
4  cg18973817  rs2476601     0   A   G   1  114471878  114377568      0

   distance                                         seq_before  \
0     53335  actgagaaagaaaattaacaagtaatacaagaaaaaactccataaa...
1     76011  GCGATGATGACGTGTGGTACTTTTCTGTTATGTTTTTAAATAATTT...
2     22880  ACTTGACAGGTCATCAGACACCCCATTACCCAAACGCCTAAAAGCC...
3     70178  ttttttttttgaaacggagtctgactctgtcgaccaggctggaatg...
4     94310  taatttttattttaggttcaggggtatgtgtgccagtttgttatgt...

                                           seq_after  seq_len
0  actgagaaagaaaattaacaagtaatacaagaaaaaactccataaa...   200001
1  GCGATGATGACGTGTGGTACTTTTCTGTTATGTTTTTAAATAATTT...   200001
2  ACTTGACAGGTCATCAGACACCCCATTACCCAAACGCCTAAAAGCC...   200001
3  ttttttttttgaaacggagtctgactctgtcgaccaggctggaatg...   200001
4  taatttttattttaggttcaggggtatgtgtgccagtttgttatgt...   200001
'''

# post-process
file_path = 'data/seq_mapping/'
output_path = 'data/seq_mapping_post/'

for snp in snp_list:
    for model in model_size.keys():
        data = pd.read_pickle(file_path + snp + '_' + model + '.dataset')
        data_check = pd.read_pickle(file_path + snp + '_' + model + '.dataset')
        max_range = model_size[model]
        for i in range(len(data_check)):
            seq_len = data_check['seq_len'][i]
            cpg = data_check['CpG'][i]
            snp = data_check['SNP'][i]
            if(seq_len != max_range * 2 + 1):
                data = data.drop(data[(data['CpG']==cpg)&(data['SNP']==snp)].index)
        data = data.reset_index(drop=True)
        data.to_pickle(output_path + snp + '_' + model + '.dataset')