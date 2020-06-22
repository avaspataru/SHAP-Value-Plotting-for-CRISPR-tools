# Shapley Value Plotting for CRISPR tools
MSc Thesis Project

This project is producing SHAP summary plots for the features of some CRISPR sgRNA design tools when ran on different datasets.

## Main files 
  **datasets** : scripts for extracting the guide sequences (in tool specific format) from the original dataset files
  
  **src/shapleyvals.py** : runs the SHAP analysis and produces a plot for the specified tool on the specified dataset
  
  **src/tooldata.py** : interface class for the tool-specific data file in order to be ran by shapleyvals.py
  
  **src/tool-model** : contains the necessary files for running the specific tool
  
## Running instructions 

```
python shapleyvals.py --tool ToolName --data DatasetName
```

Where the tool name is any of [tuscan-classification, tuscan-regression, sgRNAScorer2, wu-crispr, ssc, chop-chop-xu, chop-chop-doench, chop-chop-moreno] and dataset name is any of [xu, doench]. Case insensitive. 
  
  
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
  
  The scripts take guides of length 20 and score them assuming PAM NGG.
