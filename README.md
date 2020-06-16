# Shapley Value Plotting for CRISPR tools
MSc Thesis Project

This project is producing SHAP summary plots for the features of some CRISPR sgRNA design tools [TUSCAN, sgRNAScorer2, Wu-CRISPR, SSC] on when ran on different datasets [Xu, Doench].

## Main files 
  ./datasets : scripts for extracting the guide sequences (in tool specific format) from the original dataset files
  ./src/shapleyvals.py --tool TOOL --data DATASET : runs the SHAP analysis and produces a plot for the specified tool on the specified dataset.
  ./src/tooldata.py : interface class for the tool-specific data file in order to be ran by shapleyvals.py
  ./src/\*-model : contains the necessary files for running the specific tool
  
## To do
 - add CHOP-CHOP tool 
 - load Doench dataset in main shapleyvals
 - add other resources

## Other resources 
  - SHAP values
  - TUSCAN paper and code
  - sgRNAScorer paper and code
  - WU-CRISPR paper and code
  - SSC paper and code 
  
 ## Observations 
 For each tool, the code parts which construct the model and/or score the gRNA have been extracted and adapted to fit to the tooldata interface. For adding any other tool, its code only needs to be put into the specified format by the interface.
  
