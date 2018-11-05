# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectus T2
# Description : Performs group clssification between the four vaccine groups
# Cite        : TBD
# ******************************************************************************

rm(list = ls())

source('funcsToImport.R')

dir_res = paste('results_predClassMulti_biophysical/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

# glmnet parameters
alphas = 1  # elastic net parameter can vary between 0 (for Ridge) to 1 (for LASSO)
cvFolds = 8  # number of folds for cross-validation
repeatRun = 100 # number of repetitions of cross-validation
grpType = "ungrouped" # grouping type
intc = TRUE # intercept to the linear model?
weights_bal = TRUE # balance classes?

log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

cat('\nGlmnet Multinomial Classification parameters','\n',sep="",file=log_file,append=T)
cat('Alpha range : ',paste(alphas,'',sep=','),'\n',file=log_file,append=T)
cat('Balancing classes',weights_bal,'\n',file=log_file,append=T)
cat('Intercept for logistic regression',intc,'\n',file=log_file,append=T)
cat('Number of folds : ',cvFolds,'\n',file=log_file,append=T)
cat('Number of repeated evaluations : ',repeatRun,'\n',file=log_file,append=T)

# -------------------------------------------
# Sec 02: Data
# -------------------------------------------

luminex = read.csv('in/luminex_filtered.csv', header=TRUE, row.names=1)

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
# Sec 04: Glmnet classification
# -------------------------------------------

# Create labels for classification
classes = as.matrix(subjects[,'group_all'])
feats = luminex

dir_class = paste(dir_res,"class/",sep="")
dir.create(dir_class)

numFeat = ncol(feats)

# -------------
# Scale features and set NAs to 0
# -------------

feats = scale(feats)
na_idx = which(is.na(feats),arr.ind=TRUE)
if(length(na_idx)!=0){
  
  feats[na_idx] = 0
  
}

label = classes

weights = rep(1,length(label))

class_model = glmnetMultiClass(feats,label,weights,numFeat,grpType,intc,alphas,cvFolds,repeatRun)

# -------------------------------------------
# Sec 05: Visualize prediction performance
# -------------------------------------------

pdf(paste(dir_class,'best_model.pdf',sep=""))
plot(class_model$final_fit,main=paste('alpha: ',class_model$best_alpha,'\n',sep=""))
dev.off()

# -------------
# Plot log-odds
# -------------

pred_prob = class_model$best_model$preval
pred_class = apply(pred_prob,1,which.max)
pred_class = pred_class + 1

label_tform = numeric(length(pred_class))
for(chooseIdx in c(2,3,4,5)){
  
  label_tform[label==chooseIdx] = chooseIdx-1
  
}

log_odds = numeric(length(pred_class))
for(sampleIdx in 1:length(pred_class)){
  
  other_prob = max(pred_prob[sampleIdx,-label_tform[sampleIdx]])
  log_odds[sampleIdx] = log(pred_prob[sampleIdx,label_tform[sampleIdx]]/other_prob)
  
}

confMat = confusionMatrix(as.factor(pred_class),as.factor(label))
plotConfusionMulti(confMat,label,dir_class)

#Visualize performance
df_odds = as.data.frame(cbind(as.vector(log_odds),label,pred_class))
colnames(df_odds) = c('odds','label','label_pred')
df_odds$label = as.factor(df_odds$label)
df_odds$label_pred = as.factor(df_odds$label_pred)

pdf(paste(dir_class,'box_odds_gg.pdf',sep=""))
p1 = ggplot(df_odds,aes(x=label,y=odds)) + geom_boxplot(width=0.5,notch = F,coef=1.58,outlier.shape = NA,size=1,colour="black") + geom_point(position = position_jitter(w=0.1),size=4,aes(colour=label_pred)) + ylab('Log Odds') + scale_x_discrete(name="",labels=c('Empty','LTA1','IL-12','IL-12+LTA1')) + theme(plot.title = element_text(size=20), axis.line = element_line(colour = "black",size=1.5), axis.title.y=element_text(size=15), axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=20,colour='black'), axis.ticks=element_blank(),panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='bottom') + scale_colour_manual(values=group_colors) + scale_fill_manual(values=group_colors)
p1 = p1 + geom_hline(yintercept=0,colour='black',size=0.78,linetype='dashed',alpha=0.7)
print(p1)
dev.off()

# -------------
# Plot coefficients
# -------------

coeff_min = class_model$best_model$coeff
coeff_min_nz_idx = which(rowSums(abs(coeff_min))!=0)

pdf(paste(dir_class,'coeffs_min.pdf',sep=""))
par(mar=c(18,4,2,0.5))
barplot(t(coeff_min[coeff_min_nz_idx,]),main=c(" Predictor Coefficients"),beside=TRUE,las=2,cex.names=1.2,col=group_colors)
dev.off()

#
cat('Rept : Mean of MSE',class_model$repeat_mse_min,'(',class_model$repeat_mse_min_sd,')\n')
cat('Perm : Mean of MSE',class_model$permut_mse_min,'(',class_model$permut_mse_min_sd,')\n')

# -------------------------------------------
# Sec 07: Comparing actual and permuted models' performance
# -------------------------------------------

robust_test = as.data.frame(cbind(1-class_model$cv_repeat[,'min'],1-class_model$cv_permut[,'min']))
colnames(robust_test) = c('Luminex','Permuted')

df_test = as.data.frame(cbind(c(rep(1,repeatRun),rep(2,repeatRun)),as.vector(as.matrix(robust_test))))
colnames(df_test) = c('label','Acc')
df_test$label = as.factor(df_test$label)

diff_test = wilcox.test(robust_test[,1],robust_test[,2],alternative="two.sided")
eff_test = cliff.delta(robust_test[,1],robust_test[,2])
eff_interp = as.character(eff_test$magnitude)

cdf_perm = ecdf(robust_test[,2])
actual_avg = mean(robust_test[,1])
pval_actual = 1-cdf_perm(actual_avg)

pdf(paste(dir_class,'robust_test_v.pdf',sep=""))
p = ggplot(df_test,aes(x=label,y=Acc)) + geom_violin(size=1,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Actual','Permuted')) + scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=15,colour='black') ,axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=12,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='bottom') + scale_fill_manual(values=c('darkorange3',"#CD660066")) + xlab("") + ylab('Balanced Accuracy\n')
p = p + geom_hline(yintercept=0.25,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(robust_test$Luminex,na.rm=T),colour='darkorange3',size=1.2,linetype='dashed',alpha=0.7)
#p = p + geom_hline(yintercept=median(robust_test$Permuted,na.rm=T),colour="#CD660066",size=1.2,linetype='dashed',alpha=0.7)
# p = p + annotate("segment",x=1,xend=2,y=1,yend=1,size=2)
# p = p + annotate("text",x=1.5,y=0.97,size=6,label=paste('P : ',format(diff_test$p.value,digits=3,scientific=T)))

p = p + annotate("segment",x=1,xend=2,y=0.85,yend=0.85,size=3,colour=effect_colors[eff_interp])
p = p + annotate("rect",xmin=c(1.45), xmax=c(1.55), ymin=c(0.83) , ymax=c(0.87), color=effect_colors[eff_interp], fill=effect_colors[eff_interp])
p = p + annotate("point",shape=24,y=0.85,x=1.5,size=15,color="black",fill=effect_colors[eff_interp])
p = p + annotate("text",x=1.5,y=0.94,size=6,label=paste('P : ',format(pval_actual,digits=3,scientific=T),sep=""))

print(p)
dev.off()

# print to log file
cat('Repeated Balanced Accuracy',class_model$repeat_mse_min,'(',class_model$repeat_mse_min_sd,')\n',file=log_file,append=T)
cat('Permuted Balanced Accuracy',class_model$permut_mse_min,'(',class_model$permut_mse_min_sd,')\n',file=log_file,append=T)

cat('\n\n',rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)