from standardizeFeatures import *

nucleotides = [ 'A', 'C', 'G', 'T']
dinucleotides = [ 'AA', 'AC', 'AG', 'AT',
                'CA', 'CC', 'CG', 'CT',
                'GA', 'GC', 'GG', 'GT',
                'TA', 'TC', 'TG', 'TT']

class CoefficientModel:

    scoring = {}

    # Loads the appropriate scoring method
    # In standardized format of feature names (see standardizeFeatures.py)
    def setScoring(scoringMethod):
        if scoringMethod == 'DOENCH_2014':
            scoring = DOENCH_2014_STANDARD
        elif scoringMethod == 'XU_2015':
            scoring = XU_2015_STANDARD
        elif scoringMethod == 'MORENO_MATEOS_2015':
            scoring = MORENO_MATEOS_2015_STANDARD
        else:
            print("Scoring method not supported " + scoringMethod)
            quit()

    def predict(self,data):
        pass

class ChopChopData(ToolData):

    def loadTrainingSet(self):
        pass

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


    def loadModel(self, scoringMethod):
        model = CoefficientModel()
        model.setScoring(scoringMethod)
        return model

    def getFeatures(self,seq):
        pass
