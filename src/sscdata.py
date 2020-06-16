
import numpy as np
import pandas as pd

from tooldata import *

class MatrixModel:

    matrix = [[]]

    def loadMatrix(self):
        matrixFile = open('./ssc-model/human_CRISPRi_20bp.matrix','r')
        lines = matrixFile.readlines()
        lines.pop(0) #intercept (at 0)
        lines.pop(0) #header
        matrixFile.close()

        matrix = []
        row = 1
        for line in lines:
            coefs = line.strip().split('\t')
            coefs = [float(coef) for coef in coefs]
            matrix.append(coefs)
            row = row + 1

        self.matrix = matrix

    def predict(self, data):
        nucleotides = ['A','C','G','T']
        results = []
        for d in data:
            score = 0
            #one data point is a list of features (0-80)
            for f in range(0,20):
                for n in range(0,4):
                    #one feature is if the nucleotide is there
                    if(d[f*4 + n] == 1): #this nucleotide is present in position f
                        score = score + self.matrix[f][n]
            results.append(score)

        return np.array(results)

class SSCData(ToolData):

    def loadTrainingSet(self): #training set is Xu
        xuFile = open('../datasets/Xu-2015_Is-Efficient.csv')
        GUIDE_INDEX = 4
        GUIDE_START_POS = 10
        GUIDE_LEN = 20

        sequences = []

        # read each line - file contains header and blank footer - then break apart by comma
        for line in [x.split(',') for x in xuFile.read().split('\n')[1:-1]]:
            guideSeq = line[GUIDE_INDEX][GUIDE_START_POS:(GUIDE_START_POS + GUIDE_LEN)]
            sequences.append(guideSeq)
        xuFile.close()

        feature_set = []
        for g in sequences:
            features = self.getFeatures(g)
            feature_set.append(features)

        df = pd.DataFrame(feature_set)
        return df

    def loadFeatureNames(self):
        names = []
        nucleotides = ['A','C','G','T']
        for i in range(0,20): #features are just what is at each position (flags)
            for n in nucleotides:
                names.append(str(i) + ":" + n)
        return names


    def loadModel(self):
        model = MatrixModel()
        model.loadMatrix()
        return model

    def getFeatures(self,seq):
        features = []
        nucleotides = ['A','C','G','T']
        for s in seq:
            for n in nucleotides:
                if(s == n):
                    features.append(1)
                else:
                    features.append(0)
        return features
