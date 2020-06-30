import pickle
import argparse
import shap
import sys
import re
import matplotlib.pyplot as plt
import numpy as np

from utils import getToolObject

nucleotides = ['A','C','G','T']

#gets the tool name from the file name given as argument
def getToolName(fileName):
    toolName = fileName.split('-')[1]
    return toolName

#gets the dataset name from the fileName given as argument
def getDatasetName(fileName):
    dataName = fileName.split('-')[2]
    return dataName

#gets the feature names for the tool which was used to produce the pickle file
def getFeatureNames(fileName):

    #get the tool corresponding to the file name
    toolName = getToolName(fileName)
    tool = getToolObject(toolName) #from utils.py

    #get the actual feature names from the tool object
    featureNames = tool.loadFeatureNames()
    return featureNames

#get the pickle file from arguments
parser = argparse.ArgumentParser(description="specify the file to unpicke for the SHAP values")
parser.add_argument('--file', help='path to pickle file', required=True)
args = parser.parse_args()

fileName = args.file

#load the SHAP values and the data
p = open('../results/'+fileName,"rb")
print("Loading data...")
shapValues = pickle.load(p)
dataset = pickle.load(p)
p.close()

#feature names for the tool the pickle was created with
featureNames = getFeatureNames( fileName )

#calculate the average shap value across all datapoints for a feature
print("Calculating averages...")
avgShapVals = [0] * len(featureNames) #initialise

#shapValues = dataset x features
for d in range(0,len(dataset)):
    for f in range(0,len(featureNames)):
        avgShapVals[f] = avgShapVals[f] + shapValues[d][f] #add the value from this datapoint

#average out the sum
avgShapVals = [ v / len(dataset) for v in avgShapVals]


#sort the features by the absolute value of the avg Shap value
print("Getting positional values...")
featureImp = sorted( zip(avgShapVals,featureNames) , key = lambda t: abs(t[0])) #zip with the feature names and sort
#featureImp = featureImp[-20:] # take the top 20 highest values <=> most important 20 features
featureImp = list(reversed(featureImp)) #order from most important to least important

#create dictionary from list of feature importances
#where averages[nucleotide][pos] = the shapley value of the feature pos:nucleotide
averages = np.zeros((4,20))
matches = 0 #number of positional features found
for f in featureImp:
    featureValue = f[0] #the shapley value of this feature
    featureName = f[1] #the name of this feature

    #matching the feature name to something containing position-nucleotide [only one nucleotide!]
    #(position) (maybe :) (nucleotide) (anything but another nucleotide)
    m = re.match(r"(?P<pos>((\d)(\d)*))(:)?(?P<nuc>[ACGT])([^ACGT](.*))?$", featureName)

    if m == None: #maybe the position and nucleotide are in reverse order
        print("rev")
        #(something + not a nucleotide or start) (nucleotide) (maybe :) (pos) (not digit + something or end)
        m = re.match(r"^((.*)[^ACGT])?(?P<nuc>[ACGT])(:)?(?P<pos>((\d)(\d)*))([^0-19](.*))?$", featureName)

    if m == None: #not a positional feature
        continue
    matches = matches + 1
    pos = int( (m.group('pos')) )
    nuc = nucleotides.index(m.group('nuc')) # A:0, C:1, G:2, T:3

    averages[nuc][pos] = featureValue

assert matches == 80 #can only have 80 positional features

print("Plotting...")

width = 0.35 # the width of the bars
positions = range(0,20)
pA = plt.bar(positions, averages[0], width)
pC = plt.bar(positions, averages[1], width)
pG = plt.bar(positions, averages[2], width)
pT = plt.bar(positions, averages[3], width)

plt.xticks(positions)
plt.plot(range(-1,21), [0]*22, color='black', linewidth=0.5) #plot a line at 0

plt.ylabel("SHAP values")
plt.title("SHAP values for positional features from " + getToolName(fileName).upper() + " ran on the " + getDatasetName(fileName).upper() + " dataset")
plt.legend((pA[0], pC[0], pG[0], pT[0]), ('A', 'C', 'G', 'T'))

plt.show()
