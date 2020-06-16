from tooldata import *
import os
import subprocess, shlex
import pandas as pd
from numpy import random, array, zeros, empty
import numpy as np
import joblib

cwd = os.getcwd()

class StackingModel:

    def svm_scale(self, feature, ranges, lower, upper):
    	feature_scaled = zeros((feature.shape[0], feature.shape[1]))

    	for i in range(0, feature.shape[1]-1):
    		feat = feature[:, i]
    		min_val = ranges[i, 0]
    		max_val = ranges[i, 1]
    		if ((max_val - min_val) != 0):
    			feat_std = (feat - min_val) / (max_val - min_val)
    			feature_scaled[:, i] = feat_std * (upper - lower) + lower
    		else:
    			feature_scaled[:, i] = zeros((feature.shape[0],))

    	return feature_scaled

    def predict(self, features_test): #stacking predict from WuCrispr
        ntest = features_test.shape[0]
        second_test = zeros((ntest, 2))
        svm_test_kf = empty((5, ntest))
        boost_test_kf = empty((5, ntest))

        train_range = joblib.load( cwd + '/wu-crispr-model/models/train_range.joblib')

        features_test_scale = self.svm_scale(features_test, train_range, 0, 1)

        for i in range(0,5):
            filename1 = cwd + '/wu-crispr-model/models/SVM_model{}.joblib'.format(i)
            svm_classifier = joblib.load(filename1)
            svm_test_kf[i,:] = svm_classifier.predict_proba(features_test_scale)[:, 1]

        svm_test = svm_test_kf.mean(axis = 0)
        second_test[:, 0] = svm_test

        for j in range(0,5):
            filename2 = cwd + '/wu-crispr-model/models/XGBoost_model{}.joblib'.format(j)
            boost_classifier = joblib.load(filename2)
            boost_test_kf[j, :] = boost_classifier.predict_proba(features_test)[:, 1]

        boost_test = boost_test_kf.mean(axis = 0)
        second_test[:, 1] = boost_test

        second_filename = cwd + '/wu-crispr-model/models/second_logistic_model.joblib'
        second_layer_classifier = joblib.load(second_filename)
        label_pred = second_layer_classifier.predict(second_test)
        prob_pred = second_layer_classifier.predict_proba(second_test)[:, 1]

        return prob_pred

class WuCrisprData(ToolData):
    """A class containing the structure for a tool model data"""

    def loadTrainingSet(self):
        #the training set is the doench dataset
        doenchFile = open('../datasets/Doench-Wu-Crispr.txt') #has guide + 3 nucleotides in the tail
        guideTuples = doenchFile.readlines()
        doenchFile.close()
        guides = [ l.strip().split(' ')[0] for l in guideTuples]

        print("Loaded training set of " + str(len(guides)) + " samples.")
        print("Computing features for the training set [not already available]")
        feature_set = []
        cnt = 0
        for g in guides:
            features = self.getFeatures(g)
            feature_set.append(features)
            if cnt % 100 == 0:
                print("-- Calculated for "+ str(cnt) + " samples.")
            cnt = cnt + 1
        df = pd.DataFrame(feature_set)
        return df

    def loadFeatureNames(self):
        file = open ('./wu-crispr-model/names')
        names = file.readlines()
        names = [n.strip() for n in names]
        file.close()
        return names

    def loadModel(self):
        s = StackingModel()
        return s

    def getFeatures(self,seq):

        #call the perl script on this sequence
        command = "perl ./wu-crispr-model/callPerl.pl " + seq
        args =  shlex.split(command)
        result = subprocess.run(args, stdout=subprocess.PIPE)

        #parse output from stdout
        features = [float(i) for i in result.stdout.decode('utf-8').split(',')]
        return features
