# random forest for training enformer-classification model
import argparse
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

model_size_list = ['small','large']

final_path = '../../datasets/ten_fold/enformer/'
output_path = '../../model/pred_results_tenfold_enformer/'

for model_size in model_size_list:
    
    # load & merge datasets
    for f in range(10):
        fold = str(f+1)

        test_data = pd.read_pickle(final_path + model_size + '_split' + fold + '.dataset')
        merge_data = pd.DataFrame()
        for i in range(10):
            if(str(i+1)!= fold):
                data = pd.read_pickle(final_path + model_size + '_split' + str(i+1) + '.dataset')
                merge_data = pd.concat([merge_data,data])
        train_data = merge_data.reset_index(drop=True)

        # training model
        feature_list = []
        labels = np.array(train_data['label'].astype("int"))
        for i in range(train_data.shape[0]):
            sample_feature = []
            sample_feature += train_data['result_after'][i].flatten().tolist()
            feature_list.append(sample_feature)
        features = np.array(feature_list)

        X_train = features
        Y_train = labels
        print(X_train.shape)
        print(Y_train.shape)

        clf = RandomForestClassifier()
        clf.fit(X_train,Y_train)

        # prediction output
        feature_list = []
        labels = np.array(test_data['label'].astype("int"))
        for i in range(test_data.shape[0]):
            sample_feature = []
            sample_feature += test_data['result_after'][i].flatten().tolist()
            feature_list.append(sample_feature)
        features = np.array(feature_list)

        X_test = features
        Y_test = labels
        print(X_test.shape)
        print(Y_test.shape)

        y_score = clf.predict(X_test)                                   
        y_score_pro = clf.predict_proba(X_test) # (.., 2)

        # save results
        np.save(output_path + model_size + '_split' + fold + '_label.npy',Y_test)
        np.save(output_path + model_size + '_split' + fold + '_score.npy',y_score)
        np.save(output_path + model_size + '_split' + fold + '_score_pro.npy',y_score_pro)