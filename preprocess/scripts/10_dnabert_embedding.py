# fetch DNABERT embedding: 500bp -> [1,768]
# https://github.com/jerryji1993/DNABERT, k-mer=6
import os
import torch
from transformers import BertModel, BertConfig, DNATokenizer
import numpy as np

dir_to_pretrained_model = "/data/eqtl/dnabert/6-new-12w-0/"

config = BertConfig.from_pretrained('/home/liuzhe/dnabert/DNABERT/src/transformers/dnabert-config/bert-config-6/config.json')
tokenizer = DNATokenizer.from_pretrained('/home/liuzhe/dnabert/DNABERT/src/transformers/dnabert-config/bert-config-6/vocab.txt')
model = BertModel.from_pretrained(dir_to_pretrained_model, config=config)

'''
# test code:

sequence = "AATCTAATCTAGTCTAGCCTAGCA"
model_input = tokenizer.encode_plus(sequence, add_special_tokens=True, max_length=512)["input_ids"]
model_input = torch.tensor(model_input, dtype=torch.long)
model_input = model_input.unsqueeze(0)   # to generate a fake batch with batch size one

output = model(model_input)
new = output[1].data.cpu().numpy()
np.save('seq2_embedding.npy',new)
'''

model_size = {'2001':'small','20001':'middle','200001':'large','2000001':'huge'}
model_cutting = {'small':1,'middle':19,'large':199,'huge':1_999}

def split_seq(sequence):
    seq_list = []
    length = len(sequence)
    model = model_size[str(int(length))]
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

def dnabert_embedding(seq_list):
    embedding_list = []
    for sequence in seq_list:
        model_input = tokenizer.encode_plus(sequence, add_special_tokens=True, max_length=512)["input_ids"]
        model_input = torch.tensor(model_input, dtype=torch.long)
        model_input = model_input.unsqueeze(0)
        output = model(model_input)
        new = output[1].data.cpu().numpy()
        embedding_list.append(new)
    embedding = np.array(embedding_list[0])
    for item in embedding_list[1:]:
        embedding = np.concatenate([embedding,item],axis=0)
    return embedding

file_path = '/data/eqtl/methven/sequence_datasets/'
output_path = '/data/eqtl/methven/dnabert_embedding/'

filenames = os.listdir(file_path)
for file in filenames:
    file_clean = file.split('.')[0]
    with open(file_path + file,'r') as f1:
        sequence = f1.readline().strip()
        sequence_list = split_seq(sequence)
        dnabert = dnabert_embedding(sequence_list)
        np.save(output_path + file_clean + '.npy', dnabert)
            