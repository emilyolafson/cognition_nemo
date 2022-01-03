this repository contains code to clean the dataset from Boes et al from the Iowa Stroke Regsitry. This dataest is highly annotated and many variables are correlated.

cleandata.ipynb replaces miscellaneous characters with NaNs and produces a list of the # of NaNs per feature & per subject.
featureselect.m identifies the optimal number of features (cognitive tests) to predict (that have full coverage across all subjects) while keeping as many subjects as possible. 
PCA of these final features demonstrated that two cognitive scores represent the largest variance across the population. these variables will be predicted from structural & functional disconnectivity.

