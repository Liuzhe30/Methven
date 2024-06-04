# generate middle output of EMO layers

import argparse
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)

import os
import copy
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input

from model.methven import *
from src.dataGenerator import dataGenerator

def generate_middle_output(model_size, layer_name): 

    train_data = pd.read_pickle('datasets/final/dnabert/' + model_size + '_train.dataset')
    
    # check model size
    if(model_size == 'small'):
        batch_size = 1024
        model = build_methven_small()
    elif(model_size == 'large'):
        batch_size = 512
        model = build_methven_large()

    trainGenerator = dataGenerator(train_data, batch_size, model_size)

    model.load_weights('model/weights/' + model_size + '_trained_weights.tf').expect_partial()

    input_features = trainGenerator.generate_validation()[0]
    npy_input = input_features[0]
    npy_input2 = input_features[1]
    npy_input3 = np.expand_dims(input_features[2],axis=-1)

    npy_input = np.sum(npy_input, axis=1)
    npy_input2 = np.sum(npy_input2, axis=1)
    npy_input3 = np.sum(npy_input3, axis=1)

    layer_output = tf.keras.models.Model(inputs=model.input,outputs=model.get_layer(layer_name).output)
    npy_out = layer_output.predict(trainGenerator.generate_validation()[0], batch_size=batch_size)  
    #print(npy_out.shape)

    label = trainGenerator.generate_validation()[1]
    
    # save prediction output
    # npy_input = np.sum(npy_input, axis=1) # for large and huge model
    np.save('model/middle_output/' + model_size + '/' + 'input.npy', npy_input)
    np.save('model/middle_output/' + model_size + '/' + 'input2.npy', npy_input2)
    np.save('model/middle_output/' + model_size + '/' + 'input3.npy', npy_input3)
    np.save('model/middle_output/' + model_size + '/' + 'output.npy', npy_out)
    np.save('model/middle_output/' + model_size + '/' + 'label.npy', label)


if __name__=='__main__':

    generate_middle_output('small','dense_3')
    #generate_middle_output('large','dense_3')