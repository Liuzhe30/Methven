# class for data generation
import pandas as pd
import numpy as np
import copy
pd.set_option('display.max_columns', None)

class dataGenerator():

    def __init__(self, dataset, batch_size, model_size):
        self.dataset = dataset
        self.batch_size = batch_size
        self.model_size = model_size

        # 2-d label: one-hot       
        self.label_dict = {'0':[1,0], # no effect
                            '1':[0,1], # effect
                            }

        self.gene_dict = {'A':[0,0,0,1], 'T':[0,0,1,0], 'C':[0,1,0,0], 'G':[1,0,0,0], 
             'a':[0,0,0,1], 't':[0,0,1,0], 'c':[0,1,0,0], 'g':[1,0,0,0],
             'N':[0,0,0,0], 'n':[0,0,0,0],
             'P':[0,0,0,0]} # padding

    def generate_batch(self):

        while 1:
            i = 0
            while i < (len(self.dataset) - self.batch_size):
                before_seq, after_seq, dataY_batch = [], [], []
                for j in range(i, i + self.batch_size):

                    # generate label
                    label = str(self.dataset['label'].values[j])
                    dataY_batch.append(self.label_dict[label])

                    # generate input_51 
                    seq_list = []
                    value = self.dataset['seq_before'].values[j]
                    for strr in value:
                        seq_list.append(self.gene_dict[strr])
                    before_seq.append(seq_list)

                    seq_list = []
                    value = self.dataset['seq_after'].values[j] 
                    for strr in value:
                        seq_list.append(self.gene_dict[strr])
                    after_seq.append(seq_list)

                input_before = np.array(before_seq)
                input_after = np.array(after_seq)
                y = np.array(dataY_batch)

                i += self.batch_size
                
                yield ([input_before, input_after], y)

    def generate_validation(self):
        
        before_seq, after_seq, dataY_batch = [], [], []
        for j in range(len(self.dataset)):

            # generate label
            label = str(self.dataset['label'].values[j])
            dataY_batch.append(self.label_dict[label])

            # generate input_51 
            seq_list = []
            value = self.dataset['seq_before'].values[j]
            for strr in value:
                seq_list.append(self.gene_dict[strr])
            before_seq.append(seq_list)

            seq_list = []
            value = self.dataset['seq_after'].values[j] 
            for strr in value:
                seq_list.append(self.gene_dict[strr])
            after_seq.append(seq_list)

        input_before = np.array(before_seq)
        input_after = np.array(after_seq)
        y = np.array(dataY_batch)

        return ([input_before, input_after], y)