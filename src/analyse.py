import pickle
import argparse
import shap
import sys

def getFeatureNames(toolName):
    if toolName in ['tuscan-regression','tuscan-classification']:
        sys.path.insert(0,'../src/tuscan-model')
        from tuscandata import TuscanData
        tool = TuscanData()
        tool.setRegressionFlag(toolName=='tuscan-regression') #if regression
    elif toolName == 'sgrnascorer2':
        sys.path.insert(0, '../src/sgRNAScorer2-model')
        from sgRNAScorer2data import SgRNAScorer2Data
        tool = SgRNAScorer2Data()
    elif toolName == 'wu-crispr':
        sys.path.insert(0,'../src/wu-crispr-model')
        from wucrisprdata import WuCrisprData
        tool = WuCrisprData()
    elif toolName == 'ssc':
        sys.path.insert(0, '../src/ssc-model')
        from sscdata import SSCData
        tool = SSCData()
    elif toolName in ['chop-chop-xu', 'chop-chop-doench', 'chop-chop-moreno']:
        #ASSUME Xu scoring method = TO DO CHANGE
        sys.path.insert(0, '../src/chop-chop-model')
        from chopchopdata import ChopChopData
        tool = ChopChopData()
        tool.setScoring(toolName)
    else:
        print("Tool Name not valid")
        quit()


    featureNames = tool.loadFeatureNames()
    print(featureNames)
    quit()

#get the pickle file from arguments
parser = argparse.ArgumentParser(description="specify the file to unpicke for the SHAP values")
parser.add_argument('--file', help='path to pickle file', required=True)
args = parser.parse_args()

fileName = args.file

#load the SHAP values and the data
p = open('../results/'+fileName,"rb")
print("Loading data...")
shap_values = pickle.load(p)
dataset = pickle.load(p)
p.close()

#feature Names for the tool the pickle was created with
feature_names = getFeatureNames( fileName.split('-')[1] )
