# training

import argparse
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)

import os
import copy
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import optimizers
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input
from sklearn.utils import shuffle

from model.methven import *
from src.dataGenerator import dataGenerator

# for gpu training
# os.environ["CUDA_VISIBLE_DEVICES"] = '0' # may needed when DEVICE:0 is occupied, otherwise where will be an error about positional embedding
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print(e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model_size', default='small')
    parser.add_argument('--epochs', default=200, type=int)
    parser.add_argument('--lr', default=0.001, type=float,
                        help="Initial learning rate")
    parser.add_argument('--lr_decay', default=0.05, type=float,
                        help="The value multiplied by lr at each epoch. Set a larger value for larger epochs") 
    parser.add_argument('--save_dir', default='model/weights/')
    parser.add_argument('--debug', action='store_true',
                        help="Save weights by TensorBoard")
    args = parser.parse_args()
    model_size = args.model_size

    path = args.save_dir + model_size + '/'
    if(not os.path.exists(path)):
        os.makedirs(path)

    # load datasets
    train_data = pd.read_pickle('datasets/final/dnabert/' + model_size + '_train.dataset')
    valid_data = pd.read_pickle('datasets/final/dnabert/' + model_size + '_valid.dataset')

    # check model size
    if(model_size == 'small'):
        batch_size = 1024
        model = build_methven_small()
    elif(model_size == 'large'):
        batch_size = 512
        model = build_methven_large()

    model.summary()

    save_dir = args.save_dir
    
    # callbacks
    log = tf.keras.callbacks.CSVLogger(save_dir + model_size + '_log.csv')
    checkpoint = tf.keras.callbacks.ModelCheckpoint(save_dir + args.model_size + '_trained_weights.tf', monitor='val_acc', mode='max', #val_categorical_accuracy val_acc
                                       save_best_only=True, save_weights_only=True, verbose=1)        
    earlystop = tf.keras.callbacks.EarlyStopping(patience=10, min_delta=1e-3)

    # Train the model and save it
    model.compile(loss='binary_crossentropy', 
              optimizer='adam',
              metrics=['mae', 'acc'])
    
    trainGenerator = dataGenerator(train_data, batch_size, model_size)
    validGenerator = dataGenerator(valid_data, batch_size, model_size)
    
    history = model.fit(trainGenerator.generate_batch(), # Tf2 new feature
          steps_per_epoch = len(train_data)/batch_size,
          epochs = args.epochs, verbose=1,
          validation_data = validGenerator.generate_batch(),
          validation_steps = len(valid_data)/batch_size,
          callbacks = [log, checkpoint, earlystop],
          shuffle = True,
          workers = 1).history

    print('Trained model saved to \'%s/trained_model.tf\'' % save_dir)