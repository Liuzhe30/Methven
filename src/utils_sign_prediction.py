import sys 
sys.path.append("..") 
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input

from model.methven import *
from src.dataGenerator import dataGenerator

def get_sign_prediction_result(download_path):
    
    # 1 fetch dataset
    test_data = pd.read_pickle('temp.dataset')
    input_variant = test_data['input_variant'][0]
    cpg_position = test_data['cpg_position'][0]
    
    chr_str = input_variant.split('_')[0]
    position = int(input_variant.split('_')[1])
    before_mutation = input_variant.split('_')[2]
    after_mutation = input_variant.split('_')[3]
    distance = np.abs(int(cpg_position) - int(position))

    model_size = get_model_size(distance)
    model_max_range = {'small':10_000,'large':100_000}
    max_range = model_max_range[model_size]

    # 2 prediction
    # check model size
    if(model_size == 'small'):
        batch_size = 1
        model = build_methven_small()
    elif(model_size == 'large'):
        batch_size = 1
        model = build_methven_large()
    #model.summary()
    testGenerator = dataGenerator(test_data, batch_size, model_size)
    model.load_weights(download_path + model_size + '/weight/' + model_size + '_trained_weights.tf').expect_partial()
    results = model.predict(testGenerator.generate_validation()[0], batch_size=batch_size)
    results_dict = {}
    results_dict['score'] = results[0]
    if(np.argmax(results[0]) == 1):
        results_dict['label'] = 'Up-regulation'
    else:
        results_dict['label'] = 'Down-regulation'
    return results_dict['label']

def get_model_size(distance):
    tss_distance = np.abs(distance)
    if(tss_distance < 10000):
        model_size = 'small'
    elif(tss_distance >= 10000):
        model_size = 'large'
    return model_size