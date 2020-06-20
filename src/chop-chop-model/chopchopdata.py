from standardizeFeatures import *
import numpy as np
import pandas as pd

nucleotides = [ 'A', 'C', 'G', 'T']
dinucleotides = [ 'AA', 'AC', 'AG', 'AT',
                'CA', 'CC', 'CG', 'CT',
                'GA', 'GC', 'GG', 'GT',
                'TA', 'TC', 'TG', 'TT']

class CoefficientModel:

    scoring = {}

    # Loads the appropriate scoring method
    # In standardized format of feature names (see standardizeFeatures.py)
    def setScoring(self,scoringMethod):
        if scoringMethod == 'DOENCH_2014':
            self.scoring = getDoenchScoring()
        elif scoringMethod == 'XU_2015':
            self.scoring = getXuScoring()
        elif scoringMethod == 'MORENO_MATEOS_2015':
            self.scoring = getMorenoScoring()
        else:
            print("Scoring method not supported " + scoringMethod)
            quit()

    #For a given dataframe of data with calculated features
    #Return the scores based on the current scoring method
    def predict(self,data):
        results = []

        for d in data:
            score = 0

            #Intercept of the scoring method
            if "Intercept" in self.scoring.keys():
                score = self.scoring["Intercept"]

            #GC content features (0,1) = (low, high)
            if d[0] != 0 and "gc_low" in self.scoring.keys(): #gc_low
                gc = d[0]
                score = score + (abs(gc-10) * self.scoring["gc_low"])
            if d[1] != 0 and "gc_high" in self.scoring.keys(): #gc_high
                gc = d[1]
                score = score + ((gc-10) * self.scoring["gc_high"])

            #N in position i
            featureIndex = 2
            for i in range(0,20):
                for n in nucleotides:
                    if d[featureIndex] == 1: #yes, n is in position i
                        if n+str(i) in self.scoring.keys(): #otherwise score is 0
                            score = score + self.scoring[n+str(i)]
                    featureIndex = featureIndex + 1

            #NN in position i
            for i in range(0,19): #the guide (without last position)
                for dinuc in dinucleotides:
                    if d[featureIndex] == 1: #yes, dinuc at position i
                        if dinuc+str(i) in self.scoring.keys(): #otherwise score is 0
                            score = score + self.scoring[dinuc+str(i)]
                    featureIndex = featureIndex + 1

            #PAM nucleotides
            for i in range(0,3):
                for n in nucleotides:
                    if d[featureIndex] == 1: #yes, n in pam at pos i
                        if "PAM"+n+str(i) in self.scoring.keys(): #otherwise score is 0
                            score = score + self.scoring["PAM" + n + str(i)]
                    featureIndex = featureIndex + 1

            #PAM dinucleotides
            for i in range(0,2):
                for dinuc in dinucleotides:
                    if d[featureIndex] == 1: #yes, dinuc in pam at i
                        if "PAM"+dinuc+str(i) in self.scoring.keys(): #otherwise score is 0
                            score = score + self.scoring["PAM" + dinuc + str(i)]
                    featureIndex = featureIndex + 1

            results.append(score)

        return np.array(results)



class ChopChopData():

    #the training set is the Xu dataet for the Xu scoring method
    def loadTrainingSet(self):
        xuFile = open('../datasets/Xu-2015_Is-Efficient.csv')
        GUIDE_INDEX = 4
        GUIDE_START_POS = 10
        GUIDE_LEN = 23

        sequences = []

        # read each line - file contains header and blank footer - then break apart by comma
        for line in [x.split(',') for x in xuFile.read().split('\n')[1:-1]]:
            guideSeq = line[GUIDE_INDEX][GUIDE_START_POS:(GUIDE_START_POS + GUIDE_LEN)]
            sequences.append(guideSeq)
        xuFile.close()

        #compute the features for the training set sequences
        feature_set = []
        for g in sequences:
            features = self.getFeatures(g)
            feature_set.append(features)

        #put the features in a dataframe
        df = pd.DataFrame(feature_set)
        return df


    def loadFeatureNames(self):
        featureNames = []

        featureNames.append("gc_low") # <10
        featureNames.append("gc_high") # >10

        #nucleotide at position
        for pos in range(0,20): #the guide
            for n in nucleotides:
                featureName = n + str(pos)
                featureNames.append(featureName)

        ##dinucleotide at position
        for pos in range(0,19): #the guide (without last position)
            for dinuc in dinucleotides:
                featureName = dinuc + str(pos)
                featureNames.append(featureName)

        #PAM nucleotides
        for pamPos in range(0,3):
            for n in nucleotides:
                featureName = "PAM" + n + str(pamPos)
                featureNames.append(featureName)

        #PAM dinucleotides
        for pamPos in range(0,2):
            for dinuc in dinucleotides:
                featureName = "PAM" + dinuc + str(pamPos)
                featureNames.append(featureName)

        #since we are not evaluating any sequences with upstream dimers, no need to get those features.

        return featureNames

    def loadModel(self):
        model = CoefficientModel()
        model.setScoring("XU_2015") #TODO: put all
        return model

    def getFeatures(self,seq):
        features = []
        #separate
        pam = seq[:-3]
        guide = seq[0:20]

        #GC content features
        gc = guide.count('G') + guide.count('C')
        if gc < 10:
            features.append(gc) #gc_low
            features.append(0) #gc not high
        elif gc > 10:
            features.append(0) #gc not low
            features.append(gc) #gc high
        else:
            features.append(0) #gc not low
            features.append(0) #gc not high

        #nucleotide at position
        for pos in range(0,20): #the guide
            for n in nucleotides:
                if guide[pos] == n:
                    features.append(1)
                else:
                    features.append(0)

        ##dinucleotide at position
        for pos in range(0,19): #the guide (without last position)
            for dinuc in dinucleotides:
                if guide[pos:pos+2] == dinuc:
                    features.append(1)
                else:
                    features.append(0)

        #PAM nucleotides
        for pamPos in range(0,3):
            for n in nucleotides:
                if pam[pamPos] == n:
                    features.append(1)
                else:
                    features.append(0)

        #PAM dinucleotides
        for pamPos in range(0,2):
            for dinuc in dinucleotides:
                if pam[pamPos:pamPos+2] == dinuc:
                    features.append(1)
                else:
                    features.append(0)
        return features
