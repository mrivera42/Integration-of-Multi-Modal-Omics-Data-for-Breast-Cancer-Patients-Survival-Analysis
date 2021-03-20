# Integration-of-Multi-Modal-Omics-Data-for-Breast-Cancer-Patients-Survival-Analysis

Michael-Alexander Rivera and Shreenu Sivakumar

## Instructions:
Each of our 3 models (miRNA, mRNA, and miRNA + mRNA) has a main.m file associated with it’s name. For a chosen model, to execute all steps of our workflow, including preprocessing, feature selection, and training, and survival analysis, run the main file.

### Main Scripts:
* __main_mirna.m__ = main script for miRNA only model
* __main_rnaseq.m__ = main script for mRNA only model 
* __main_integrated.m__ = main script for miRNA + mRNA integrated model

### Supporting Functions and Scripts:
* __preprocessing_mirna.m__ = performs preprocessing for miRNA data (no inputs needed) 
* __preprocessing_rnaseq.m__ = performs preprocessing for mRNA data (no inputs needed) 
* __preprocessing_survival.m__ = performs preprocessing for survival data (no inputs needed) 
* __featureselection_mirna.m__ = performs feature selection on preprocessed data (either modality) 
* __cox_mirna.m__ = calculates risk scores for patients in testing data (either modality) 
* __kaplanmeier.m__ = uses risk scores to split patients into high and low risk groups for visualization 
* __concordanceIndex.m__ = uses survival time and risk scores to calculate concordance index 
* __stepFunction.m__ = supporting function for concordance index
* __tables.m__ = script that create heat maps and plots to organize performance metric results 

### Data Files used:
* __mirna.tsv__ = tsv file of miRNA data 
* __rnaseq-aaa.tsv__ = tsv file of mRNA 
* __data survival.tsv__ = tsv file for survival data
