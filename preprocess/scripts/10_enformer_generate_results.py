# generate enformer results for comparison (enformer-cpu)
# dimensionality reduction using TSNE-PCA 
import os
import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow_hub as hub
from sklearn.manifold import TSNE

# https://github.com/google-deepmind/deepmind-research/blob/master/enformer/README.md, one hot encoded in order 'ACGT'
eyes = np.eye(4)
gene_dict = {'A':eyes[0], 'C':eyes[1], 'G':eyes[2], 'T':eyes[3], 'N':np.zeros(4),
             'a':eyes[0], 'c':eyes[1], 'g':eyes[2], 't':eyes[3], 'n':np.zeros(4)
             }

enformer = hub.load("https://www.kaggle.com/models/deepmind/enformer/frameworks/TensorFlow2/variations/enformer/versions/1").model
print('test load successful!')
maxlen = 393216
fasta_path = '../../datasets/raw/reference_genome_hg19/'

model_size = {'small':10_000,'large':100_000}

def mutation_center_seq(chr, cpg_pos, snp_pos, after_mutation):
    chr_str = chr
    with open(fasta_path + chr_str + '.fasta') as fa:
        line = fa.readline()        
        range_seq = int(maxlen/2)       
        sequence_before = line[cpg_pos - range_seq:cpg_pos + range_seq]
        sequence_after = line[cpg_pos - range_seq:snp_pos - 1] + after_mutation + line[snp_pos:cpg_pos + range_seq]
    return sequence_before, sequence_after

def fetch_enformer_results(sequence):
    seq_list = []
    for strr in sequence:
        seq_list.append(gene_dict[strr])
    seq_array = np.array(seq_list)
    tensor = tf.convert_to_tensor(seq_array, tf.float32)
    tensor = tf.expand_dims(tensor, axis=0)
    stack = enformer.predict_on_batch(tensor)['human']
    stack = np.array(stack)
    new_stack = stack.reshape((stack.shape[1],stack.shape[2]))
    # generate pca results
    tsne = TSNE(n_components=3, init='pca', random_state=0)
    result = tsne.fit_transform(new_stack)

    return result

for i in range(10):
    for model_size in model_size.keys():
        fold = str(i+1)
        train_all = pd.read_pickle('../../datasets/ten_fold/origin/' + model_size + '_split' + fold + '.dataset')
        train_df = pd.DataFrame(columns=['CpG', 'SNP', 'label', 'result_before', 'result_after'])

        for i in range(train_all.shape[0]):
            chr = train_all['CHR'].values[i]
            cpg_pos = train_all['CpG_POS'].values[i]
            snp_pos = train_all['SNP_POS'].values[i]
            after_mutation = train_all['Alt'].values[i]
            sequence_before, sequence_after = mutation_center_seq(chr, cpg_pos, snp_pos, after_mutation)
            if(len(sequence_before) != maxlen or len(sequence_after) != maxlen):
                continue
            result_before = fetch_enformer_results(sequence_before)
            result_after = fetch_enformer_results(sequence_after)
            print(result_before.shape)
            train_df = train_df._append([{'CpG': train_all['CpG'][i], 'SNP': train_all['SNP'][i], 'label': train_all['label'][i], 'result_before': result_before, 
                                            'result_after': result_after}], ignore_index=True)
        train_df.to_pickle('../../datasets/ten_fold/enformer/' + model_size + '_split' + fold + '.dataset')



