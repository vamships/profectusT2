# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectus T2
# Description : Imports the necessary packages, functions, and defines global variables
# Cite        : TBD
# ******************************************************************************

library(caret)
library(corrplot)
library(e1071)
library(effsize)
library(ggplot2)
library(glmnet)
library(gplots)
library(polycor)
library(survcomp)
library(survival)

# -----------------------------------------------------
source('funcs/convertKMtoEvents.R')
source('funcs/coxSurvivalFinalModel.R')
source('funcs/coxSurvivalLookBack.R')
source('funcs/coxSurvivalWithBackSearch.R')
source('funcs/createColumnColors.R')
source('funcs/createSubjectColors.R')
source('funcs/doCoxBackSearch.R')
source('funcs/doFullSurvivalAnalysis.R')
source('funcs/extractProbabilityFromKM.R')
source('funcs/glmnetBiClass.R')
source('funcs/glmnetMultiClass.R')
source('funcs/heatmap4.R')
source('funcs/plotConfusion.R')
source('funcs/takeOffOneFeat.R')
source('funcs/univariatePolyserialFilter.R')
# -----------------------------------------------------
# Define colors

# group type
group_id = 2:5
group_colors = c('skyblue2','dodgerblue2','firebrick','purple4')
names(group_colors) = group_id

# colors based on how many challenges they survived
challenge_colors = c(colorRampPalette(c('honeydew2','grey34'))(10), 'lawngreen');
names(challenge_colors) = c(1:10, 'UI')

# the two groups of interest : IL12 and Others
km_colors = c('blue4','firebrick')
names(km_colors) = c(1,2)

# effect size
effect_colors = c('grey83','grey63','grey43','grey13')
names(effect_colors) = c('negligible','small','medium','large')

# antigen and reagent names
antigen_names = c('SIV.gp160.Ertl.Ag','SIVcpz.EK505.gp120','SIVmac.gp140.Novartis.P173.Ag','SIVsmH4.PR55.Gag','SIVmac239.BK28.PR55','SIVmac239.gp120','SIVmac239.gp130','SIVmac239.gp140','SIVsmE543.3.cv2a','SIVsmE543.3.cv2c','SIVsmE543.gp140','SIVmac1A11.gp140','3352.FLSC','CCG7V.rhFLSC','SHIV162P3.gp120','SIVmac239.Pol..Ecoli.','mac239.rhFLSC')
reagent_names = c('IgG.Low',"aRhIgG.PE.low",'IgG.High',"aRhIgG.PE.high",'R2A.1','R2A.2','R2A.3','R2A.4','R3A.1','R3A.3','FcgRIIa','FcgRIIIa','FcgRIIIb','C1q','MBL','R2B.1')

# antigen colors
antigen_colors = NULL
antigen_colors['SIV.gp160.Ertl.Ag'] = "snow3"
antigen_colors['SIVcpz.EK505.gp120'] = "blue2"
antigen_colors['SIVsmH4.PR55.Gag'] = "cadetblue2"
antigen_colors['SIVmac239.BK28.PR55'] = "cadetblue2"
antigen_colors['SIVmac239.gp120'] = "greenyellow"
antigen_colors['SIVmac239.gp130'] = "greenyellow"
antigen_colors['SIVmac239.gp140'] = "greenyellow"
antigen_colors['SIVsmE543.3.cv2a'] = "darkgoldenrod4"
antigen_colors['SIVsmE543.3.cv2c'] = "darkgoldenrod4"
antigen_colors['SIVsmE543.gp140'] = "lightgoldenrod3"
antigen_colors['SIVmac1A11.gp140'] = "darkviolet"
antigen_colors['3352.FLSC'] = "gold"
antigen_colors['CCG7V.rhFLSC'] = "goldenrod2"
antigen_colors['mac239.rhFLSC'] = "khaki"
antigen_colors['SHIV162P3.gp120'] = "navajowhite"
antigen_colors['SIVmac239.Pol..Ecoli.'] = "lightslateblue"
antigen_colors['SIVmac.gp140.Novartis.P173.Ag'] = "brown2"

# reagent colors
reagent_colors = NULL
reagent_colors['IgG.Low'] = "peru"
reagent_colors['aRhIgG.PE.low'] = "peru"
reagent_colors['aRhIgG.PE.high'] = "saddlebrown"
reagent_colors['IgG.High'] = "saddlebrown"
reagent_colors['R2A.1'] = "springgreen4"
reagent_colors['R2A.2'] = "springgreen4"
reagent_colors['R2A.3'] = "springgreen4"
reagent_colors['R2A.4'] = "springgreen4"
reagent_colors['R3A.1'] = "wheat4"
reagent_colors['R3A.3'] = "wheat4"
reagent_colors['FcgRIIa'] = "lightseagreen"
reagent_colors['FcgRIIIa'] = "khaki4"
reagent_colors['FcgRIIIb'] = "violetred4"
reagent_colors['C1q'] = "chocolate2"
reagent_colors['MBL'] = "seashell3"
reagent_colors['R2B.1'] = "thistle1"

# reagent_list_1 = c('IgG.Low','IgG.High','R2A.1','R3A.1','R3A.3','FcgRIIa','FcgRIIIa','FcgRIIIb','C1q','MBL')
# reagent_list_2 = c('aRhIgG.PE.low','aRhIgG.PE.high','R2A.2','R2A.3','R2A.4','R2B.1','R3A.1','R3A.3','C1q','MBL')
# 
# antigen_list_1 = c('SIV.gp160.Ertl.Ag','SIVcpz.EK505.gp120','SIVmac.gp140.Novartis.P173.Ag','SIVsmH4.PR55.Gag','SIVmac239.BK28.PR55','SIVmac239.gp120','SIVmac239.gp130','SIVmac239.gp140','SIVsmE543.3.cv2a','SIVsmE543.3.cv2c','SIVsmE543.gp140','SIVmac1A11.gp140')
# antigen_list_2 = c('3352.FLSC','CCG7V.rhFLSC','SHIV162P3.gp120','SIVcpz.EK505.gp120','SIVmac239.Pol..Ecoli.','SIVmac239.gp120','SIVmac239.gp140','SIVsmH4.p55.Gag','mac239.rhFLSC')