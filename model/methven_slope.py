# model structure 

import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)

import os
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input

def attention_3d_block(inputs, nb_time_steps):
    a = layers.Permute((2, 1))(inputs)
    a = layers.Dense(nb_time_steps, activation='relu', name = 'attention_dense')(a)
    a_probs = layers.Permute((2, 1), name='attention_vec')(a)
    #output_attention_mul = merge([inputs, a_probs], name='attention_mul', mode='mul')
    output_attention_mul = layers.multiply([inputs, a_probs], name='attention_mul')
    return output_attention_mul

def build_methven_small():

    # hyper-paramaters
    maxlen = 41

    ####### inputs (from data generator)
    input1 = Input(shape=(maxlen, 768), name = 'input_before') # seq-before-mutation
    input2 = Input(shape=(maxlen, 768), name = 'input_after') # seq-after-mutation
    input3 = Input(shape=(maxlen,), name = 'input_atac') # atac

    ####### merge inputs of the same scale
    new_input3 = layers.Reshape((maxlen,1))(input3)
    input_between = layers.concatenate([input1, input2, new_input3], axis=-1)
    print('input_between.get_shape()', input_between.get_shape()) # (None, 5, 1537)
    
    ####### between-branch
    attention_mul = attention_3d_block(input_between,maxlen)
    attention_mul = layers.BatchNormalization()(attention_mul)
    attention_mul = layers.Dropout(0.3)(attention_mul)

    gru_out = layers.Bidirectional(layers.GRU(64, activation='tanh', recurrent_activation='sigmoid', use_bias=True, return_sequences=True))(input_between)
    gru_out = layers.Bidirectional(layers.GRU(64, activation='tanh', recurrent_activation='sigmoid', use_bias=True, return_sequences=True))(gru_out)
    gru_out = layers.BatchNormalization()(gru_out)
    #print('gru_out.get_shape()', gru_out.get_shape()) # (None, maxlen, 64)

    fc = layers.Dense(128, activation='relu',kernel_initializer='random_uniform', bias_initializer='zeros')(gru_out)  
    fc = layers.BatchNormalization()(fc)
    fc = layers.Dense(32, activation='relu')(fc)     
    fc = layers.Dense(16, activation='relu')(fc)
    fc = layers.Flatten()(fc)
    fc = layers.Dense(16)(fc)
    output = layers.Dense(1)(fc)
    
    model = Model(inputs=[input1, input2, input3], outputs=output)
    return model


def build_methven_large():

    # hyper-paramaters
    maxlen = 400

    ####### inputs (from data generator)
    input1 = Input(shape=(maxlen, 768), name = 'input_before') # seq-before-mutation
    input2 = Input(shape=(maxlen, 768), name = 'input_after') # seq-after-mutation
    input3 = Input(shape=(maxlen,), name = 'input_atac') # atac

    ####### merge inputs of the same scale
    new_input3 = layers.Reshape((maxlen,1))(input3)
    input_between = layers.concatenate([input1, input2, new_input3], axis=-1)
    print('input_between.get_shape()', input_between.get_shape()) # (None, 5, 1537)
    
    ####### between-branch
    attention_mul = attention_3d_block(input_between,maxlen)
    attention_mul = layers.BatchNormalization()(attention_mul)
    attention_mul = layers.Dropout(0.3)(attention_mul)

    gru_out = layers.Bidirectional(layers.GRU(64, activation='tanh', recurrent_activation='sigmoid', use_bias=True, return_sequences=True))(input_between)
    gru_out = layers.Bidirectional(layers.GRU(64, activation='tanh', recurrent_activation='sigmoid', use_bias=True, return_sequences=True))(gru_out)
    gru_out = layers.BatchNormalization()(gru_out)
    #print('gru_out.get_shape()', gru_out.get_shape()) # (None, maxlen, 64)

    fc = layers.Dense(128, activation='relu',kernel_initializer='random_uniform', bias_initializer='zeros')(gru_out)  
    fc = layers.BatchNormalization()(fc)
    fc = layers.Dense(32, activation='relu')(fc)     
    fc = layers.Dense(16, activation='relu')(fc)
    fc = layers.Flatten()(fc)
    fc = layers.Dense(16)(fc)
    output = layers.Dense(1)(fc)
    
    model = Model(inputs=[input1, input2, input3], outputs=output)
    return model