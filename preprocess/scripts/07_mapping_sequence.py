# mapping DNA sequences, centered on CpG
import pandas as pd
pd.set_option('display.max_columns', None)

fasta_path = '../../datasets/raw/reference_genome_hg38/'
file_path = '../../datasets/chr_downsampled/'
output_path = '../../datasets/seq_mapping/'

'''
        CpG         SNP    Beta Ref Alt CHR   CpG_POS   SNP_POS distance label                                                                                                                                                                                                                                | 6/2000 [00:01<07:03,  4.71it/s] 
0  cg10676837  rs10158705 -0.6795   G   A   1  19452577  19424086    28491     1
1  cg17081867  rs10158705 -0.4874   G   A   1  19441601  19424086    17515     1
2  cg08558153  rs10158705 -0.2824   G   A   1  19446335  19424086    22249     1
3  cg25159064  rs10158705 -0.2571   G   A   1  19390988  19424086    33098     1
4  cg01832549  rs10158705  0.2954   G   A   1  19448494  19424086    24408     1
'''

model_size = {'small':1_000,'middle':10_000,'large':100_000,'huge':1_000_000}
for i in range(22):
    chr = 'chr' + str(i+1)
    for model in model_size.keys():
        data = pd.read_pickle(file_path + chr + '_' + model + '.dataset')
        data['seq_before'] = ''
        data['seq_after'] = ''
        data['seq_len'] = 0
        max_range = model_size[model]
        for i in range(len(data)):
            after_mutation = data['Alt'].values[i]
            cpg_pos = int(data['CpG_POS'].values[i])
            snp_pos = int(data['SNP_POS'].values[i])

            with open(fasta_path + chr + '.fasta') as fa:
                line = fa.readline()
                seq_before = line[cpg_pos - max_range - 1:cpg_pos + max_range]
                seq_len = len(seq_before)
                seq_after = line[cpg_pos - max_range - 1:snp_pos - 1] + after_mutation + line[snp_pos:cpg_pos + max_range] # replace A1

            data.loc[i, 'seq_before'] = seq_before
            data.loc[i, 'seq_after'] = seq_after
            data.loc[i, 'seq_len'] = seq_len
            #print(data.head())
        data.to_pickle(output_path + chr + '_' + model + '.dataset')

'''
          CpG         SNP    Beta Ref Alt CHR    CpG_POS    SNP_POS distance  \
0  cg27024417  rs72666775  0.0000   G   A   1   32180211   32179560      651
1  cg22302152   rs4660159  0.0000   T   C   1   40097901   40098028      127
2  cg16008418   rs6692840  0.5934   A   G   1   64186457   64186821      364
3  cg11650704   rs1127314 -0.3977   G   A   1  154584098  154583790      308
4  cg10220895   rs1188722  0.0000   T   C   1  228276435  228276575      140

  label                                         seq_before  \
0     0  cctttggaagaaggacattgaggcccagagagagaacagaacgtcc...
1     0  CTCTCTCTTTCCTTCCCCTTCTCTCTTTCCAGTGAAGGGGACGGCT...
2     1  ggaatttttcagctccattagtatcttttgggactaccgttttaaa...
3     1  GCTGGCTCAGCACGTTCCCAGATGAAGGTCTGTCATCTAAAGGAGA...
4     0  AGTGCTACTAGTGGGGATCACGTCTAGGATTTCGAGGCACAAGGGC...

                                           seq_after  seq_len
0  cctttggaagaaggacattgaggcccagagagagaacagaacgtcc...     2001
1  CTCTCTCTTTCCTTCCCCTTCTCTCTTTCCAGTGAAGGGGACGGCT...     2001
2  ggaatttttcagctccattagtatcttttgggactaccgttttaaa...     2001
3  GCTGGCTCAGCACGTTCCCAGATGAAGGTCTGTCATCTAAAGGAGA...     2001
4  AGTGCTACTAGTGGGGATCACGTCTAGGATTTCGAGGCACAAGGGC...     2001
'''