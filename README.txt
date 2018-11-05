# ******************************************************************************
# Author      : Srivamshi Pittala (srivamshi.pittala.gr@dartmouth.edu)
# Advisor     : Prof. Chris Bailey-Kellogg (cbk@cs.dartmouth.edu) & Prof. Margaret E. Ackerman
# Project     : Profectus T2
# Description : Help file for running the scripts to analyze luminex and functional measurements
#				associated with the Profectus HIV vaccine trial T2
# Cite        : TBD
# ******************************************************************************

#***********************
Input Files (in the directory in/)
1. luminex_filtered.csv      ||
2. luminex_titer.csv         || Biophysical measurements
3. luminex_titerAdjusted.csv ||
4. func_response.csv 	: Functional and T-cell measurements
5. subjects.csv		: group and survival information
#***********************

#***********************
For survival analysis, use the following three scripts
#***********************

#-----------------------
1. coxPredRisk_biophysical.R
Takes about 5 min to complete
#-----------------------

===> This script performs the survival analysis using the biophysical measurements (luminex). 
		(1)	First, the features are pre-filtered using the polyserial correlation coefficient.
		(2) Then repeated cross-validation is performed on the features returned from pre-filtering. In this 	cross-validation, greedy backward elimination is used to select features for each fold using the training set. The most-frequent features appearing in the repeated cross-validation are used to build a final model, the results from which are used for the figures.
		(3) The features in the final model are used to discover other features that could be equally predictive of risk but were not considered since they are correlated with the final feature set.
		(4) Permutation tests are performed by shuffling the rows of the outcome labels, but keeping the rows of feature matrix the same. The permuted data are sent through the same pipeline as were the actual data (i.e. pre-filtering, repeated cross-validation, and final model evaluation). This is repeated multiple times, independently permuting the outcome labels every time.
		(5) Robustness is estimated by comparing the C-indices from using actual features to those of using permuted features.

#-----------------------
2. coxPredRisk_titer.R
Takes about 5 min to complete
#-----------------------

===> This script performs the same survival analysis as above, but using the biophysical measurements corresponding to antibody titer.

#-----------------------
3. coxPredRisk_titerAdjusted.R
Takes about 5 min to complete
#-----------------------

===> This script performs the same survival analysis as above, but using the biophysical measurements that are adjusted to antibody titer.

#***********************
For binomial logistic classification, use the following three scripts
#***********************

#-----------------------
1. predClassBinary_biophysical.R
Takes about 5 min to complete
#-----------------------
===> This script performs binomial logistic classification to identify vaccine groups using the biophysical measurements
		(1) Classification is done using the lasso-regularized binomial logistic regression on the features. The best regularization parameter (lambda) is chosen to be the one with lowest classification error. This is repeated multiple times. A final model is trained and evaluated by using a fixed seed to determine folds.
		(2) Permutation tests are performed by shuffling the rows of the class labels, but keeping the rows of feature matrix the same. The permuted data are sent through the same pipeline as were the actual data. This is repeated multiple times, independently permuting the class labels every time.
		(3) Robustness is estimated by comparing the accuracies from using actual features to those of using permuted features.

#-----------------------
2. predClassBinary_titer.R
Takes about 5 min to complete
#-----------------------

===> This script performs the same classification as above, but using the biophysical measurements corresponding to antibody titer.

#-----------------------
2. predClassBinary_titerAdjusted.R
Takes about 5 min to complete
#-----------------------

===> This script performs the same classification as above, but using the biophysical measurements that are adjusted to antibody titer.

#***********************
For comparing the risk predictions with functional and T-cell responses, use the following script
#***********************

#-----------------------
1. visualize_func_relation.R
Takes about 1 min to complete
#-----------------------

===> This script reads the risk predictions from coxPredRisk_biophysical.R and compares with functional(ADCC) and T-cell(IFN-gamma) measurements.

#***********************
For multinomial logistic classification, use the following script
#***********************

#-----------------------
1. predClassMulti_biophysical.R
Takes about 5 min to complete
#-----------------------
===> This script performs multinomial logistic classification to identify vaccine groups using the biophysical measurements
		(1) Classification is done using the lasso-regularized multinomial logistic regression on the features. The best regularization parameter (lambda) is chosen to be the one with lowest classification error. This is repeated multiple times. A final model is trained and evaluated by using a fixed seed to determine folds.
		(2) Permutation tests are performed by shuffling the rows of the class labels, but keeping the rows of feature matrix the same. The permuted data are sent through the same pipeline as were the actual data. This is repeated multiple times, independently permuting the class labels every time.
		(3) Robustness is estimated by comparing the accuracies from using actual features to those of using permuted and random features.

#***********************
For generating the figures as in the manuscript, use the following two scripts
#***********************

#-----------------------
1. generate_networkPlot.R
Takes about 2 min to complete
#-----------------------
===> This script generates the network plot to visualize the co-correlates of the biophysical measurements identified by the final model. Run this after 'coxPredRisk_biophysical.R' has successfully finished without errors.
	(1) The nodes of the network are the features. Edges are used to representation correlations that are above 0.75
	(2) The output of the file is an xml file, that can be used by a network visualization software like Cytoscape to finally generate the figure.
	(3) I used Cytoscape software (version 3.5.1). The steps follow:
		(i)		Open the Cytoscape software application
		(ii)	File->Import->Network->File...
		(iii)	Choose the nodes.xml file
		(iv)	Layout->Attribute Circle Layout->group
		(v)		File->Export as Image...
		(vi)	Save file as Features.png in the same directory as nodes.xml

#-----------------------
2. generate_figures.R
#-----------------------
===> This script generates the figures as shown in the manuscript. Run this after after all the above scripts have finished successfully. The directory ‘results_figures_reference/’ can be used as a reference to what the final figure outputs should look like.

# ******************************************************************************

#-----------------------
System configuration 
#-----------------------
OS 	: Ubuntu 14.04.5 LTS
CPU : i7 8 cores @ 3.6 GHz
RAM : 24 GB

#-----------------------
Software and packages (version)
#-----------------------
01. R (3.4.3) and RStudio (1.0.143)
02. glmnet (2.0-13)
03. gplots (3.0.1)
04. ggplot2 (2.2.1)
05. effsize (0.7.1)
06. survival (2.41-3)
07. corrplot (0.84)
08. caret (6.0-77)	(Install with dependencies=TRUE)
09. survcomp (1.24.0)	(Installed via bioconductor)
10. polycor (0.7-9)
11. e1071(1.6-8)

#-----------------------
Functions used by the scripts (in the directory funcs/)
#-----------------------
01. convertKMtoEvents.R
02. coxSurvivalFinalModel.R
03. coxSurvivalLookBack.R
04. coxSurvivalWithBackSearch.R
05. createColumnColors.R
06. createSubjectColors.R
07. doCoxBackSearch.R
08. doFullSurvivalAnalysis.R
09. extractProbabilityFromKM.R
10. glmnetBiClass.R
11. glmnetMultiClass.R
12. heatmap4.R
13. plotConfusion.R
14. takeOffOneFeat.R
15. univariatePolyserialFilter.R
