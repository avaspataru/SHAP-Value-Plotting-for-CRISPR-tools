# Shapley Value Plotting for CRISPR tools
MSc Thesis Project

This project is producing SHAP summary plots for the features of some CRISPR sgRNA design tools when ran on different datasets. The main functionality is that given a tool name and a dataset name, it will run the SHAP analysis and produce a plot with the SHAP values for the top 20 most important features. The values are also saved in pickle files to be loaded more easily afterwards. 

## Setting up 
 The scripts in this project require Python 3 and the following packages (which can be installed with pip3) with their corresponding dependencies: 
 
 ```
joblib==0.15.1
matplotlib==3.2.2
pandas==1.0.4
pybedtools==0.8.1
scikit-learn==0.20.0
shap==0.35.0
```
For a fast set-up create a conda environment with the provided environment.yml file. This contains all the packages needed by the project. For instructions on how to load and use an environment, refer to: [conda env](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

## Running instructions 

**To run the complete SHAP analysis**:

```
python computeShapVals.py --tool ToolName --data DatasetName
```

Where the tool name is any of [tuscan-classification, tuscan-regression, sgRNAScorer2, wu-crispr, ssc, chop-chop-xu, chop-chop-doench, chop-chop-moreno] and dataset name is any of [xu, doench]. Case insensitive. 

The plot will appear on the screen and the SHAP values will be saved in results/SHAP-toolName-datasetName as a pickle file. This is an example output:

![alt text](https://github.com/avaspataru/Shapley-Value-Plotting-for-CRISPR-tools/blob/master/plots/SHAP-summary/Tuscan-Classification-Doench.JPG)
  
**To compute the plot from the presaved SHAP values**:

```
python plotShapFromPickle.py --file pickleFileName
```

Where the pickleFileName must be a pickle file in the results directory and contain two pickles (one with the SHAP values and one with the data). These files are produced automatically by the computeShapVals.py script. 

The plot produced will be the same as the above one, but the running time will be much faster. 

**To produce plots with the SHAP values for positional features**: 

```
python plotPositionsFromPickle.py --file pickleFileName
```

Where the pickleFileName is the same as above (from the results directory). This script will produce a plot looking at the positional features (e.g. G at position 19 in the guide) and plot the SHAP values in a bar plot.

The plot will appear on the screen and will have the SHAP values for each position in the guide. This is an example output: 

![alt text](https://github.com/avaspataru/Shapley-Value-Plotting-for-CRISPR-tools/blob/master/plots/guide-positions/SSC-xu.JPG)

**To compare positions across different methods**: 

```
python comparePositionTools.py
```

This script will produce a plot similar to the one above (for one tool, all the positional features), but by combining multiple tools. This serves to compare the positional preferences across tools. It requires the pickle files to be generated for the tools to be used.

**To create a heatmap and csv of all the average SHAP values across all the tools**: 

```
python shapToHeatMap.py
```

This script will take the average SHAP values of all positional features across all the tools and prints them in a csv (avg_shap_vals.csv) and produces a heatmap showing all the values together. It requires the pickle files for all tools to be generated. 

## Main files 
  **datasets** : scripts for extracting the guide sequences (in tool specific format) from the original dataset files
  
  **src/shapleyvals.py** : runs the SHAP analysis and produces a plot for the specified tool on the specified dataset 
  
  **src/tooldata.py** : interface class for the tool-specific data file in order to be ran by shapleyvals.py
  
  **src/tool-model** : contains the necessary files for running the specific tool
  
  **results**: contains pickle files with saved SHAP values for all the models and tools ran
  
  **plotfrompickle.py**: given a pickle file name (from within results), it will load the SHAP values and produce the SHAP summary plot
  
  **analyse.py**: given a pickle file name (from within results), it will generate a plot of SHAP values by positions in the guide
  
  **utils.py**: this contains general methods used by the other scripts

  
## Other resources 
  - [Explain your model with the shap values](https://towardsdatascience.com/explain-your-model-with-the-shap-values-bc36aac4de3d)
  - TUSCAN [paper](https://pubmed.ncbi.nlm.nih.gov/31021206/) and [code](https://github.com/BauerLab/TUSCAN)
  - sgRNAScorer [paper](https://pubmed.ncbi.nlm.nih.gov/28146356/) and [code](https://sgrnascorer.cancer.gov/)
  - WU-CRISPR [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0) and [code](https://github.com/wang-lab/WU-CRISPR)
  - SSC [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4509999/) and [code](https://sourceforge.net/projects/spacerscoringcrispr/) 
  - CHOP-CHOP [paper](https://academic.oup.com/nar/article/47/W1/W171/5491735) and [code](https://bitbucket.org/valenlab/chopchop/src/master/)
 ## Observations 
 For each tool, the code parts which construct the model and/or score the gRNA have been extracted and adapted to fit to the tooldata interface. For adding any other tool, its code only needs to be put into the specified format by the interface.
 
 Some of the tools require different packages to unpickle files (different versions of scikit-learn). There will be a warning informing you if the wrong version is ran.
  
  The scripts take guides of length 20 and score them assuming PAM NGG. All the guide positions are standardized to 0-19, where 19 is the position adjacent to the PAM.
