# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectus T2
# Description : Performs survival analysis using the titer-adjusted biophysical measurements
# Cite        : TBD
# ******************************************************************************

# Copyright (C) <2018>  <Srivamshi Pittala>

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

rm(list = ls())

source('funcsToImport.R')

dir_res = paste('results_coxPred_titerAdjusted/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

num_folds = 8 # number of folds for cross-validation
featSelectMethod = 'upf' # Univariate polyserial filter
thresh_upf = 0.9 # polyserial coefficient threshold
stop_limit = 0.3 # log-likelihood cutoff for stopping backward search
num_repeat = 100 # number of repetitions of cross-validation
top_feat_thresh = 90 # Frequency cutoff to determine most-frequent features
seed_final = 8357 # seed to determine folds for final cross-validation model
doLookBack = F # To look for correlates of correlates after the final model
lookBackCorThresh = 0.75 # correlation cutoff to determine if features are correlated
set_plots = TRUE

log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

cat('\nCox-model parameters','\n',sep="",file=log_file,append=T)
cat('Number of folds : ',num_folds,'\n',file=log_file,append=T)
cat('Feature Selection Method : ',featSelectMethod,'\n',file=log_file,append=T)
cat('Threshold for upf : ',thresh_upf,'\n',file=log_file,append=T)
cat('Stop Limit for Backward Elimiation : ',stop_limit,'\n',file=log_file,append=T)
cat('Number of repeated evaluations : ',num_repeat,'\n',file=log_file,append=T)
cat('Threshold for frequency of features : ',top_feat_thresh,'\n',file=log_file,append=T)
cat('Seed for final model : ',seed_final,'\n',file=log_file,append=T)
cat('Activate lookBack feature : ',doLookBack,'\n',file=log_file,append=T)
cat('Correlation threshold for looking back : ',lookBackCorThresh,'\n',file=log_file,append=T)

# -------------------------------------------
# Sec 02: Data
# -------------------------------------------

luminex = read.csv(paste('in/luminex_titerAdjusted.csv',sep=""), header=TRUE, row.names=1)

subjects = read.csv('in/subjects.csv', header=TRUE, row.names=1)

featNames = colnames(luminex)
lcolors = createColumnColors(featNames,reagent_names,reagent_colors,antigen_names,antigen_colors)

scolors = createSubjectColors(subjects,group_colors,km_colors,challenge_colors)

run_color = rep("chocolate",ncol(luminex))
run_color[grep('_r2',colnames(luminex))] = "cornflowerblue"

lcolors = cbind(lcolors,run_color)

# -------------------------------------------
# Sec 03: Legend
# -------------------------------------------

pdf(paste(dir_res,'legend.pdf',sep=""))
plot.new()
legend("topleft",legend=names(reagent_colors),fill=reagent_colors,cex=0.66)
legend("top",legend=names(antigen_colors),fill=antigen_colors,cex=0.66)
legend("right",legend=c('Empty','LTA1','IL12','LTA1+IL12'),fill=group_colors,cex=0.66)
legend("left",legend=names(challenge_colors),fill=challenge_colors,cex=0.66,horiz = T)
legend("bottomleft",legend=c('+ve','-ve'),fill=c('khaki4','darkgoldenrod2'),cex=1.2)
legend("bottomright",legend=c('Predicted','Observed'),lty=c(1,2),cex=1.2)
dev.off()

# -------------------------------------------
# Sec 04: Feature Prefiltering
# -------------------------------------------

if(featSelectMethod=='upf'){
  
  dir_featSel = paste(dir_res,'upf/',sep='')
  dir.create(dir_featSel)
  
  selected_feat_idx = univariatePolyserialFilter(thresh_upf,luminex,subjects,scolors,lcolors,plots=TRUE,dir_featSel)
  
}

cat(selected_feat_idx,sep=",",'\n')

cat('Final selected feature size - ',length(selected_feat_idx),'\n')

feats = data.frame(luminex[,selected_feat_idx])
lcolors_original = lcolors[selected_feat_idx,]

cat('\n\nFeature filtering\n',file=log_file,append=T)
cat('Feature size after filtering: ',ncol(feats),'\n',file=log_file,append=T)

# -------------------------------------------
# Sec 05: Survival Analysis
# -------------------------------------------

dir_surv = paste(dir_res,"surv/",sep="")
dir.create(dir_surv)

surv_performance = doFullSurvivalAnalysis(feats,subjects,selected_feat_idx,num_folds,stop_limit,num_repeat,top_feat_thresh,seed_final,scolors,group_colors,km_colors,lcolors_original,set_plots,dir_surv)

actual_final = surv_performance$cindex_final
actual_repeat = surv_performance$cindex_repeat

# -------------
# Lookback at features for other correlates
# -------------

if(doLookBack){
  
  cat('\n\nLookBack feature active\n',file=log_file,append=T)
  dir_lookBack = paste(dir_res,'lookback/',sep="")
  dir.create(dir_lookBack)
  
  coxSurvivalLookBack(feats,luminex,surv_performance$top_feat_idx,lookBackCorThresh,actual_final,seed_final,subjects,scolors,lcolors,group_colors,km_colors,num_folds,log_file,dir_lookBack)
  
}

# -------------------------------------------
# Sec 06: Permutation test
# -------------------------------------------

dir_perm = paste(dir_res,"perm/",sep="")
dir.create(dir_perm)

perm_repeat = matrix(NA,nrow=num_repeat,ncol=1)
colnames(perm_repeat) = 'cindex_test'

perm_log = paste(dir_perm,'perms.txt',sep="")
file.create(perm_log)

for(testIdx in 1:num_repeat){
  
  cat(rep("#",30),"\n",sep="")
  cat("Permutation test :",testIdx,"\n")
  
  dir_test = paste(dir_perm,"test_",testIdx,"/",sep="")
  dir.create(dir_test)
  
  sample_size = nrow(luminex)
  
  perm_order = sample(sample_size)
  feats_perm = data.frame(luminex[perm_order,])
  
  cat(testIdx,':',perm_order,'\n',file=perm_log,append=T)
  
  # -------------
  # Feature Prefiltering
  # -------------
  
  if(featSelectMethod=='upf'){
    
    dir_featSel = paste(dir_test,'upf/',sep='')
    dir.create(dir_featSel)
    
    selected_feat_idx = univariatePolyserialFilter(thresh_upf,feats_perm,subjects,scolors,lcolors,plots=FALSE,dir_featSel)
    
  }
  
  cat(selected_feat_idx,sep=",",'\n')
  
  cat('Final selected feature size - ',length(selected_feat_idx),'\n')
  
  feats_x = feats_perm[,selected_feat_idx]
  
  lcolors_perm = lcolors[selected_feat_idx,]
  
  # -------------
  # Survival Analysis
  # -------------
  
  model_perm = coxSurvivalWithBackSearch(feats_x,subjects,selected_feat_idx,group_colors,num_folds,stop_limit,plots=FALSE,dir_test)
  
  perm_repeat[testIdx,'cindex_test'] = model_perm$cindex_test
  
}

# -------------------------------------------
# Sec 07: Comparing actual and permuted models' performance
# -------------------------------------------

cdf_perm = ecdf(perm_repeat[,'cindex_test'])
actual_avg = mean(actual_repeat)
pval_actual = 1-cdf_perm(actual_avg)
perm_repeat = data.frame(perm_repeat)

# -------------
# Percentile of actual final model in permuted models
# -------------

pdf(paste(dir_perm,'compare_holistic.pdf',sep=""))
p = ggplot(perm_repeat,aes(x=cindex_test)) + geom_density(aes(y=..scaled..)) + scale_x_continuous(limits=c(0.3,1)) + ggtitle('Distribution of Permuted Models\n') + ylab('Density (scaled)\n') + theme(plot.title = element_text(size=25), axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=20,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=15,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),aspect.ratio=1, legend.position='bottom') + xlab('\nC-index')
p = p + geom_vline(xintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + annotate("text",x=0.45,y=0.25,size=6,label='Random',color='black')

p = p + geom_vline(xintercept=actual_avg,colour='darkorange3',size=0.78,linetype='dashed',alpha=0.7)
p = p + annotate("text",x=actual_avg+0.06,y=0.8,size=6,label=paste('Actual\n p : ',round(pval_actual,digits=2),sep=""),color='darkorange3')

p = p + geom_vline(xintercept=median(perm_repeat[,'cindex_test'],na.rm=T),colour='deeppink4',size=0.78,linetype='dashed',alpha=0.7)
p = p + annotate("text",x=median(perm_repeat[,'cindex_test'],na.rm=T),y=0.5,size=6,label='Permuted\nMedian',color='deeppink4')

print(p)
dev.off()

# -------------
# Comparing repeated cross-validation between actual and permuted
# -------------

df_test = as.data.frame(cbind(c(rep(1,num_repeat),rep(2,num_repeat)),as.vector(as.matrix(cbind(actual_repeat,perm_repeat[,'cindex_test'])))))
colnames(df_test) = c('label','C_index')
df_test$label = as.factor(df_test$label)

write.csv(df_test,file=paste(dir_perm,'robust.csv',sep=""),row.names = T)

diff_test = ks.test(actual_repeat,perm_repeat[,"cindex_test"],alternative="two.sided")

eff_test = cliff.delta(actual_repeat,perm_repeat[,"cindex_test"])
eff_interp = as.character(eff_test$magnitude)

pdf(paste(dir_perm,'robust_test.pdf',sep=""))
p = ggplot(df_test,aes(x=label,y=C_index)) + geom_violin(size=1,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Actual','Permuted')) + scale_y_continuous(limits = c(0.25,1), breaks=seq(0.25,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=12,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='bottom') + scale_fill_manual(values=c('darkorange3',"#CD660066")) + xlab("") + ylab('Concordance Index\n')
p = p + geom_hline(yintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(actual_repeat,na.rm=T),colour='darkorange3',size=1.2,linetype='dashed',alpha=0.7)
#p = p + geom_hline(yintercept=mean(perm_repeat[,'cindex_test'],na.rm=T),colour="#CD660066",size=1.2,linetype='dashed',alpha=0.7)
#p = p + geom_point(aes(x=1,y=actual_avg),shape=18,size=5,color='darkblue')
#p = p + geom_hline(yintercept=actual_avg,colour='darkblue',size=0.7,linetype='dashed',alpha=0.7)
p = p + annotate("segment",x=1,xend=2,y=0.95,yend=0.95,size=3,colour=effect_colors[eff_interp])
p = p + annotate("point",shape=24,y=0.95,x=1.5,size=15,color="black",fill=effect_colors[eff_interp])
p = p + annotate("text",x=1.5,y=0.9,size=6,label=paste('P : ',format(pval_actual,digits=3,scientific=T)))
print(p)
dev.off()

cat('\nModel Evaluation results','\n',file=log_file,append=T)
cat('Concordance Index\n',file=log_file,append=T)
cat('Repeated Cox : ',actual_final,'(Final model)',actual_avg,'(Test)','\n',file=log_file,append=T)
cat('Permuted Cox : ',mean(perm_repeat[,1],na.rm=T),'\n',file=log_file,append=T)

cat('\n\n',rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)
