# evaluation

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

def softmax(vec):
    """Compute the softmax in a numerically stable way."""
    vec = vec - np.max(vec)  # softmax(x) = softmax(x+c)
    exp_x = np.exp(vec)
    softmax_x = exp_x / np.sum(exp_x)
    return softmax_x

def evaluate(target_model, data):
	_, acc = target_model.evaluate(data)
	print("Restore model, accuracy: {:5.2f}%".format(100*acc))

def predicting(model_size): 

    test_data = pd.read_pickle('datasets/final/dnabert/' + model_size + '_test.dataset')
    
    # check model size
    if(model_size == 'small'):
        batch_size = 1024
        model = build_methven_small()
    elif(model_size == 'large'):
        batch_size = 512
        model = build_methven_large()
    model.summary()

    testGenerator = dataGenerator(test_data, batch_size, model_size)

    model.load_weights('model/weights/' + model_size + '/' + model_size + '_trained_weights.tf').expect_partial()
    results = model.predict(testGenerator.generate_validation()[0], batch_size=batch_size)
    label = testGenerator.generate_validation()[1]
    
    # save prediction output
    np.save('model/pred_results/' + model_size + '/' + 'predict.npy', results)
    np.save('model/pred_results/' + model_size + '/' + 'label.npy', label)
    

if __name__=='__main__':

    predicting('small')
    #predicting('middle')
    #predicting('large')