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
                before_seq, after_seq, atac_bet, dataY_batch = [], [], [], []
                for j in range(i, i + self.batch_size):

                    # generate label
                    label = float(self.dataset['Beta'].values[j])
                    dataY_batch.append(label)

                    # generate input sequence-between
                    value = np.array(self.dataset['dnabert_before'].values[j])[0]
                    before_seq.append(value)

                    value = np.array(self.dataset['dnabert_after'].values[j])[0]
                    after_seq.append(value)

                    # generate input atac-between
                    atac_between = list(self.dataset['atac_between'].values[j])
                    seq_list = []
                    for item in atac_between:
                        seq_list.append(item)
                    atac_bet.append(seq_list)

                input_before = np.array(before_seq)
                input_after = np.array(after_seq)
                input_atac_bet = np.array(atac_bet)
                y = np.array(dataY_batch)

                i += self.batch_size
                
                yield ([input_before, input_after, input_atac_bet], y)

    def generate_validation(self):
        
        before_seq, after_seq, atac_bet, dataY_batch = [], [], [], []
        for j in range(len(self.dataset)):

            # generate label
            label = float(self.dataset['Beta'].values[j])
            dataY_batch.append(label)

            # generate input_51 
            value = np.array(self.dataset['dnabert_before'].values[j])[0]
            before_seq.append(value)

            value = np.array(self.dataset['dnabert_after'].values[j])[0]
            after_seq.append(value)

            # generate input_bet
            atac_between = list(self.dataset['atac_between'].values[j])
            seq_list = []
            for item in atac_between:
                seq_list.append(item)
            atac_bet.append(seq_list)

        input_before = np.array(before_seq)
        input_after = np.array(after_seq)
        input_atac_bet = np.array(atac_bet)
        y = np.array(dataY_batch)

        return ([input_before, input_after, input_atac_bet], y)