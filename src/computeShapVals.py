import shap
import pandas as pd
import pickle
import argparse
import sys
import numpy as np

from utils import getToolObject
from utils import getDataset


#runs the Shapley Value analysis on the requested dataset and tool
def main(toolName,datasetName):

    tool = getToolObject(toolName)

    #load the dataset [used as data in shapley plot]
    sequences = getDataset(datasetName, toolName)
    print("Loaded the dataset " + datasetName + " of " + str(len(sequences)) + " data points.")

    #Calculate the features for the dataset
    print("Calculating features for the dataset. [Might take a while for wu-crispr]")
    feature_set = []
    cnt = 1
    for seq in sequences:
        features = tool.getFeatures(seq)
        feature_set.append(features)
        if toolName == 'wu-crispr' and cnt % 100 == 0: #inform on progress [wu-crispr takes a while]
            print("-- Calculated features for " + str(cnt))
        cnt = cnt + 1
    print("Calculated the features for this dataset.")

    #Get feature names [for printing in the shapley plot]
    feature_names = tool.loadFeatureNames()
    print("Loaded the names of all " + str(len(feature_names)) + " features.")

    #Put together features with names in one dataframe
    dataset_df = pd.DataFrame(np.array(feature_set),columns=feature_names)

    #training set (must be loaded before model)
    train_df = tool.loadTrainingSet()
    print("Loaded the training set used for tool " + toolName + ", size: " + str(train_df.shape)+ ".")

    #load model of the tool
    model = tool.loadModel()
    print("Loaded the model for tool " + toolName + ".")

    #summarize training set and subsample test set
    summary_train_df = shap.kmeans(train_df,2)
    dataset_sub_df = dataset_df # optional to speed up things can use dataset_df.sample(400)

    #compute and plot shapley values
    shap_explainer = shap.KernelExplainer(model.predict,summary_train_df)
    shap_values = shap_explainer.shap_values(dataset_sub_df)

    #save the values
    with open("../results/SHAP-"+toolName+"-"+datasetName, 'wb') as file:
        pickle.dump(shap_values,file) #the SHAP values
        pickle.dump(dataset_sub_df,file) #the data used
    print("Computed and saved SHAP values.")

    #ploatting
    shap.summary_plot(shap_values, dataset_sub_df)



#get the tool name and dataset name to start analysis

parser = argparse.ArgumentParser(description="specify the tool name and dataset name.")
parser.add_argument('--tool', help='name of the tool to be analysed', required=True)
parser.add_argument('--data', help='name of the dataset to run the tool on', required=True)
args = parser.parse_args()

toolName = args.tool.lower()
datasetName = args.data.lower()
main(toolName,datasetName)
