import sys 
sys.path.append("..") 
import torch
import pandas as pd
import numpy as np
from transformers import AutoTokenizer, AutoModel
import warnings
warnings.filterwarnings('ignore')

tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
model_dnabert = AutoModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)

def get_embedding(input_variant, cpg_position, atac_seq, genome_path, save_file):

    # 1 prepare dataset
    chr_str = input_variant.split('_')[0]
    position = int(input_variant.split('_')[1])
    before_mutation = input_variant.split('_')[2]
    after_mutation = input_variant.split('_')[3]
    distance = np.abs(int(cpg_position) - int(position))

    model_size = get_model_size(distance)
    model_max_range = {'small':10_000,'large':100_000}
    max_range = model_max_range[model_size]
    with open(genome_path + chr_str + '.fasta') as fa:
        line = fa.readline()
        seq_before = line[cpg_position - max_range - 1:cpg_position + max_range]
        seq_len = len(seq_before)
        seq_after = line[cpg_position - max_range - 1:cpg_position - 1] + after_mutation + line[position:cpg_position + max_range] # replace A1
    
    # get DNABert2 results
    sequence1_list = split_seq(seq_before, model_size)
    dnabert1 = dnabert_embedding(sequence1_list) # (x,768)
    sequence2_list = split_seq(seq_after, model_size)
    dnabert2 = dnabert_embedding(sequence2_list) # (x,768)

    # average pooling for atac-seq
    atac_pooling = atac_average_pooling(model_size, atac_seq)

    test_data = pd.DataFrame(columns=['input_variant', 'cpg_position', 'dnabert_before', 'dnabert_after', 'atac_between', 'label', 'Beta'])
    test_data = test_data._append([{'input_variant':input_variant, 'cpg_position':cpg_position, 'dnabert_before':np.expand_dims(dnabert1,axis=0), 'dnabert_after':np.expand_dims(dnabert2,axis=0), 
                                    'atac_between':atac_pooling , 'label':0, 'Beta':0}], ignore_index=True)
    
    test_data.to_pickle(save_file)

def get_model_size(distance):
    tss_distance = np.abs(distance)
    if(tss_distance < 10000):
        model_size = 'small'
    elif(tss_distance >= 10000):
        model_size = 'large'
    return model_size

def split_seq(sequence,model):

    model_cutting = {'small':19,'large':199}
    seq_list = []
    n = model_cutting[model]

    # first 250bp
    seq_list.append(sequence[0:250])
    # first n*500bp
    for i in range(n):
        seq_list.append(sequence[250+i*500:250+(i+1)*500])
    # mutation part
    seq_list.append(sequence[250+500*n:250+500*n+501])
    # second n*500bp
    for i in range(n):
        seq_list.append(sequence[250+500*n+501+i*500:250+500*n+501+(i+1)*500])
    # second 250bp
    seq_list.append(sequence[250+500*n+501+500*n:250+500*n+501+500*n+250])
    
    return seq_list

def atac_average_pooling(model_size, atac):
    model_cutting = {'small':19,'large':199}
    new_atac_list = []
    n = model_cutting[model_size]
    # first 250bp
    new_atac_list.append(np.mean(atac[0:250]))
    # first n*500bp
    for i in range(n):
        new_atac_list.append(np.mean(atac[250+i*500:250+(i+1)*500]))
    # mutation part
    new_atac_list.append(np.mean(atac[250+500*n:250+500*n+501]))
    # second n*500bp
    for i in range(n):
        new_atac_list.append(np.mean(atac[250+500*n+501+i*500:250+500*n+501+(i+1)*500]))
    # second 250bp
    new_atac_list.append(np.mean(atac[250+500*n+501+500*n:250+500*n+501+500*n+250]))
    return new_atac_list

def dnabert_embedding(seq_list):
    embedding_list = []
    for sequence in seq_list:
        inputs = tokenizer(sequence, return_tensors = 'pt')["input_ids"]
        hidden_states = model_dnabert(inputs)[0] # [1, sequence_length, 768]
        embedding_mean = torch.mean(hidden_states[0], dim=0)
        new = embedding_mean.data.cpu().numpy()
        new = new.reshape([1,new.shape[0]])
        embedding_list.append(new)
    embedding = np.array(embedding_list[0])
    for item in embedding_list[1:]:
        embedding = np.concatenate([embedding,item],axis=0)
    return embedding