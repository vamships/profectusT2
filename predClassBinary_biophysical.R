# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectus T2
# Description : Performs group clssification between IL12 and Others group
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

dir_res = paste('results_predClassBinary_biophysical/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

# glmnet parameters
alphas = 1  # elastic net parameter can vary between 0 (for Ridge) to 1 (for LASSO)
cvFolds = 8  # number of folds for cross-validation
repeatRun = 100 # number of repetitions of cross-validation
intc = TRUE # intercept to the linear model?
weights_bal = TRUE # balance classes?

log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

cat('\nGlmnet Binomial Classification parameters','\n',sep="",file=log_file,append=T)
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
classes = as.matrix((subjects[,'groupID']==2)*1)
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

cat('\n\nClassification Results ','\n',file=log_file,append=T)
cat('\nPredictions for groups classes\n',file=log_file,append=T)

label = classes
weights = rep(1,length(label))
weights[which(label==1)] = sum(label==0)/sum(label==1)

class_model = glmnetBiClass(feats,label,weights,numFeat,intc,alphas,cvFolds,repeatRun)

# -------------------------------------------
# Sec 05: Visualize prediction performance
# -------------------------------------------

pdf(paste(dir_class,'best_model.pdf',sep=""))
plot(class_model$final_fit,main=paste('alpha: ',class_model$best_alpha,'\n',sep=""))
dev.off()

# -------------
# Plot log-odds
# -------------

pred_prob = class_model$final_fit$fit.preval[,match(class_model$final_fit$lambda.min,class_model$final_fit$lambda)]
pred_class = (pred_prob>0.5)*1
confMat = confusionMatrix(as.factor(pred_class),as.factor(label))
plotConfusionBinary(confMat,label,dir_class)

df = as.data.frame(cbind(label,pred_prob,subjects[,'group_all']))
colnames(df) = c('label','lp','Group')
df$label = as.factor(df$label)
df$Group = as.factor(df$Group)
pred_table = confMat$table
class_tot = colSums(pred_table)
pdf(paste(dir_class,'box_pred.pdf',sep=""))
p = ggplot(df, aes(x=label, y=lp, color=Group)) + geom_boxplot(width=0.3,notch = F,outlier.shape = NA, na.rm=T, size=1,colour="black") + geom_point(position = position_jitter(w=0.1),size=4.5,aes(colour=Group)) + scale_x_discrete(labels=c('Others','IL-12')) + ylab('Probability of being IL-12') + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=25,colour='black'), axis.title.y = element_text(size=25,colour='black') ,axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=20,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),aspect.ratio=1, legend.position='bottom') + xlab('') + scale_colour_manual(values=group_colors)
p = p + geom_hline(yintercept=0.5,colour='black',size=1.1,linetype='dashed')
p = p + annotate("text",x=0.64,y=0.25,size=8,label=paste("frac(",pred_table[1,1],",",class_tot[1],")",sep=""),parse=T)
p = p + annotate("text",x=0.64,y=0.75,size=8,label=paste("frac(",pred_table[2,1],",",class_tot[1],")",sep=""),parse=T)
p = p + annotate("text",x=2.34,y=0.75,size=8,label=paste("frac(",pred_table[2,2],",",class_tot[2],")",sep=""),parse=T)
p = p + annotate("text",x=2.34,y=0.25,size=8,label=paste("frac(",pred_table[1,2],",",class_tot[2],")",sep=""),parse=T)
print(p)
dev.off()

# -------------
# Plot coefficients
# -------------

coeff_min = class_model$coeff_min
coeff_min_idx = class_model$coeff_min_idx
coeff_min_nz_idx = which(coeff_min!=0)

pdf(paste(dir_class,'coeffs_min.pdf',sep=""))
if(length(coeff_min_nz_idx)!=0){
  
  coeff_min = coeff_min[coeff_min_nz_idx]
  coeff_min_idx = coeff_min_idx[coeff_min_nz_idx]
  
  coeff_colors = rep(km_colors[1],length(coeff_min))
  coeff_colors[which(coeff_min>0)] = km_colors[2]
  
  par(mar=c(17,7,2,0.5))
  barplot(coeff_min,col=coeff_colors,main=c("Predictor Coefficients"),las=2,names.arg=names(coeff_min),cex.names=1.4,ylim=c(-1,1),width=0.75,ylab='',cex.axis = 2)
  mtext(expression(paste('Coefficient',sep="")), side=2, line=4, cex=3)
  
  
  
}else{
  
  plot(c(0,1),c(0,1),type='n',xaxt='n',yaxt='n',xlab='',ylab=''); text(0.5,0.5,'(empty model)')
  
}
dev.off()

coeffs_sel = coeff_min_idx[c(1,length(coeff_min))]

feats_sel = luminex[,coeffs_sel]

pdf(paste(dir_class,"biplot.pdf",sep=""))
df = data.frame(cbind(feats_sel,subjects[,'group_all']))
colnames(df) = c('f1','f2','label')
df$label = as.factor(df$label)

p = ggplot(df,aes(x=f2,y=f1,colour=label)) + geom_point(size=4,aes(colour=label)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=20,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=15,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='none') + scale_color_manual(values=group_colors) + scale_fill_manual(values=group_colors) + xlab(colnames(feats_sel)[2]) + ylab(colnames(feats_sel)[1])
print(p)
dev.off()

#
cat('Rept : Mean of Classification Error',class_model$repeat_mse_min,'(',class_model$repeat_mse_min_sd,')\n')
cat('Perm : Mean of Classification Error',class_model$permut_mse_min,'(',class_model$permut_mse_min_sd,')\n')

robust_test = as.data.frame(cbind(1-class_model$cv_repeat[,'min'],1-class_model$cv_permut[,'min']))
colnames(robust_test) = c('Actual','Permuted')

df_test = as.data.frame(cbind(c(rep(1,repeatRun),rep(2,repeatRun)),as.vector(as.matrix(robust_test))))
colnames(df_test) = c('label','Acc')
df_test$label = as.factor(df_test$label)

diff_test = wilcox.test(robust_test[,1],robust_test[,2],alternative="two.sided")
eff_test = cliff.delta(robust_test[,1],robust_test[,2])
eff_interp = as.character(eff_test$magnitude)

write.csv(df_test,file=paste(dir_class,'robust.csv',sep=""),row.names = T)

pdf(paste(dir_class,'robust_test_v.pdf',sep=""))
p = ggplot(df_test,aes(x=label,y=Acc)) + geom_violin(size=1,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Actual','Permuted')) + scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=15,colour='black') ,axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=12,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='bottom') + scale_fill_manual(values=c('darkorange3','darkslategray4')) + xlab("") + ylab('Balanced Accuracy\n')
p = p + geom_hline(yintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=median(robust_test$Actual,na.rm=T),colour='darkorange3',size=1.2,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=median(robust_test$Permuted,na.rm=T),colour='darkslategray4',size=1.2,linetype='dashed',alpha=0.7)
p = p + annotate("segment",x=1,xend=2,y=0.1,yend=0.1,size=2)
p = p + annotate("text",x=1.5,y=0.05,size=6,label=paste('P : ',format(diff_test$p.value,digits=3,scientific=T)))
p = p + annotate("text",x=1.5,y=0.15,size=6,label=paste('eff : ',format(eff_test$estimate,digits=3,scientific=T),eff_interp))
print(p)
dev.off()

# print to log file
cat('Repeated Balanced Accuracy',class_model$repeat_mse_min,'(',class_model$repeat_mse_min_sd,')\n',file=log_file,append=T)
cat('Permuted Balanced Accuracy',class_model$permut_mse_min,'(',class_model$permut_mse_min_sd,')\n',file=log_file,append=T)

cat('\n\n',rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)
