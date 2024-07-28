# cpgenie structure: https://github.com/gifford-lab/CpGenie/blob/master/cnn/seq_128x3_5_5_2f_simple.template

import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)

import os
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.constraints import MaxNorm
from tensorflow.keras import layers
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input

def build_cpgenie_small():

    # hyper-paramaters
    maxlen = 20_001
    DROPOUT = 0.3
    W_maxnorm = 3

    ####### inputs (from data generator)
    input1 = Input(shape=(maxlen, 4), name = 'input_before') # seq-before-mutation
    input2 = Input(shape=(maxlen, 4), name = 'input_after') # seq-after-mutation

    ####### merge inputs of the same scale
    input_between = layers.concatenate([input1, input2], axis=-1)
    input_between = layers.Reshape((1,input_between.get_shape()[1],input_between.get_shape()[2]))(input_between)
    print('input_between.get_shape()', input_between.get_shape()) # (None, 1, 20001, 8)

    ####### build model
    cnn = layers.Convolution2D(128, 1, 5, padding='same',activation='relu',kernel_constraint=MaxNorm(W_maxnorm))(input_between)
    cnn = layers.MaxPooling2D(pool_size=(1, 5),strides=(1,3))(cnn)
    cnn = layers.Convolution2D(256,1,5, padding='same',activation='relu',kernel_constraint=MaxNorm(W_maxnorm))(cnn)
    cnn = layers.MaxPooling2D(pool_size=(1, 5),strides=(1,3))(cnn)
    cnn = layers.Convolution2D(512,1,5, padding='same',activation='relu',kernel_constraint=MaxNorm(W_maxnorm))(cnn)
    cnn = layers.MaxPooling2D(pool_size=(1, 5),strides=(1,3))(cnn)
    print('cnn.get_shape()', cnn.get_shape()) # (None, 1, 5, 512)
    cnn = layers.Flatten()(cnn)

    fc = layers.Dense(64,activation='relu')(cnn)
    fc = layers.Dropout(DROPOUT)(fc)
    fc = layers.Dense(64,activation='relu')(fc)
    fc = layers.Dropout(DROPOUT)(fc)
    fc = layers.Dense(2)(fc)
    fc = layers.Activation('softmax')(fc)

    output = layers.Dense(2, activation = 'softmax')(fc)
    
    model = Model(inputs=[input1, input2], outputs=output)
    return model

def build_cpgenie_large():

    # hyper-paramaters
    maxlen = 200_001
    DROPOUT = 0.3
    W_maxnorm = 3

    ####### inputs (from data generator)
    input1 = Input(shape=(maxlen, 4), name = 'input_before') # seq-before-mutation
    input2 = Input(shape=(maxlen, 4), name = 'input_after') # seq-after-mutation

    ####### merge inputs of the same scale
    input_between = layers.concatenate([input1, input2], axis=-1)
    input_between = layers.Reshape((1,input_between.get_shape()[1],input_between.get_shape()[2]))(input_between)
    print('input_between.get_shape()', input_between.get_shape()) # (None, 200001, 8)

    ####### build model
    cnn = layers.Convolution2D(128, 1, 5, padding='same',activation='relu',kernel_constraint=MaxNorm(W_maxnorm))(input_between)
    cnn = layers.MaxPooling2D(pool_size=(1, 5),strides=(1,3))(cnn)
    cnn = layers.Convolution2D(256,1,5, padding='same',activation='relu',kernel_constraint=MaxNorm(W_maxnorm))(cnn)
    cnn = layers.MaxPooling2D(pool_size=(1, 5),strides=(1,3))(cnn)
    cnn = layers.Convolution2D(512,1,5, padding='same',activation='relu',kernel_constraint=MaxNorm(W_maxnorm))(cnn)
    cnn = layers.MaxPooling2D(pool_size=(1, 5),strides=(1,3))(cnn)
    cnn = layers.Flatten()(cnn)

    fc = layers.Dense(64,activation='relu')(cnn)
    fc = layers.Dropout(DROPOUT)(fc)
    fc = layers.Dense(64,activation='relu')(fc)
    fc = layers.Dropout(DROPOUT)(fc)
    fc = layers.Dense(2)(fc)
    fc = layers.Activation('softmax')(fc)

    output = layers.Dense(2, activation = 'softmax')(fc)
    
    model = Model(inputs=[input1, input2], outputs=output)
    return model