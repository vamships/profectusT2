# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectus T2
# Description : Performs survival analysis using the biophysical measurements
# Cite        : TBD
# ******************************************************************************

rm(list = ls())

source('funcsToImport.R')

dir_res = paste('results_figures/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

num_folds = 8 # number of folds for cross-validation
num_repeat = 100 # number of repetitions of cross-validation
seed_final = 8357 # seed to determine folds for final cross-validation model
set_plots = TRUE

# Survival analysis
dir_cph_biophysical = 'results_coxPred_biophysical/'
dir_cph_titer = 'results_coxPred_titer/'
dir_cph_titerAdj = 'results_coxPred_titerAdjusted/'

dir_networkPlot = 'results_networkPlot/'

# Group classification : Binary
dir_class_biophysical = 'results_predClassBinary_biophysical/'
dir_class_titer = 'results_predClassBinary_titer/'
dir_class_titerAdj = 'results_predClassBinary_titerAdjusted/'

# Group classification : Multi
dir_multi_class = 'results_predClassMulti_biophysical/'

# Tcell relation
dir_tcell_relation = 'results_func_relation/'

log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

# -------------------------------------------
# Sec 02: Data
# -------------------------------------------
luminex = read.csv('in/luminex_filtered.csv', header=TRUE, row.names=1)
luminex_adj = read.csv('in/luminex_titerAdjusted.csv', header=TRUE, row.names=1)

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
# Sec 04: Figure 2
# -------------------------------------------

dir_fig_2 = paste(dir_res,'fig_2/',sep="")
dir.create(dir_fig_2)

file.copy(paste(dir_cph_biophysical,'surv/final/km_compare.pdf',sep=""),paste(dir_fig_2,'fig_f2a.pdf',sep=""))
file.copy(paste(dir_cph_biophysical,'surv/final/CvR_test.pdf',sep=""),paste(dir_fig_2,'fig_f2b.pdf',sep=""))
file.copy(paste(dir_cph_biophysical,'surv/final/CvR_compare_test_2.pdf',sep=""),paste(dir_fig_2,'fig_f2c.pdf',sep=""))
file.copy(paste(dir_cph_biophysical,'perm/robust_test.pdf',sep=""),paste(dir_fig_2,'fig_f2d.pdf',sep=""))
file.copy(paste(dir_cph_biophysical,'surv/final/final_feat_selection_sorted.pdf',sep=""),paste(dir_fig_2,'fig_f2e.pdf',sep=""))
file.copy(paste(dir_networkPlot,'network.png',sep=""),paste(dir_fig_2,'fig_f2f.png',sep=""))

robust = read.csv(paste(dir_cph_biophysical,'perm/robust.csv',sep=""),row.names = 1)
temp = cbind(rep(NA,100),rep(2,100))
colnames(temp) = colnames(robust)

robust = rbind(robust,temp)
robust$label = as.factor(robust$label)
cindex_actual = robust[1:num_repeat,2]
cindex_perm = robust[(num_repeat+1):(num_repeat*2),2]

diff_test = wilcox.test(cindex_actual,cindex_perm,alternative="two.sided")
eff_test = cliff.delta(cindex_actual,cindex_perm)
eff_interp = as.character(eff_test$magnitude)

cdf_perm = ecdf(cindex_perm)
actual_avg = mean(cindex_actual)
pval_actual = 1-cdf_perm(actual_avg)

pdf(paste(dir_fig_2,'fig_f2d.pdf',sep=""))
p = ggplot(robust,aes(x=label,y=C_index)) + geom_violin(size=1,width=0.6,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Actual','Permuted')) + scale_y_continuous(limits = c(0.25,1), breaks=seq(0.25,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=25,colour='black') ,axis.text.x = element_text(size=25,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='bottom') + scale_fill_manual(values=c('darkorange3',"#CD660066")) + xlab("") + ylab('Concordance Index')
p = p + geom_hline(yintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(cindex_actual,na.rm=T),colour='darkorange3',size=1.2,linetype='dashed',alpha=0.7)
p = p + annotate("segment",x=1,xend=2,y=0.95,yend=0.95,size=3,colour=effect_colors[eff_interp])
p = p + annotate("point",shape=24,y=0.95,x=1.5,size=15,color="black",fill=effect_colors[eff_interp])
p = p + annotate("text",x=1.5,y=0.9,size=6,label=paste('P : ',format(pval_actual,digits=3,scientific=T),sep=""))
print(p)
dev.off()

# -------------------------------------------
# Sec 04: Figure 3
# -------------------------------------------

dir_fig_3 = paste(dir_res,'fig_3/',sep="")
dir.create(dir_fig_3)

file.copy(paste(dir_cph_titer,'surv/final/km_compare.pdf',sep=""),paste(dir_fig_3,'fig_f3a.pdf',sep=""))
file.copy(paste(dir_cph_titer,'surv/final/CvR_test.pdf',sep=""),paste(dir_fig_3,'fig_f3b.pdf',sep=""))
file.copy(paste(dir_cph_titer,'surv/final/CvR_compare_test_2.pdf',sep=""),paste(dir_fig_3,'fig_f3c.pdf',sep=""))

file.copy(paste(dir_cph_titerAdj,'surv/final/km_compare.pdf',sep=""),paste(dir_fig_3,'fig_f3d.pdf',sep=""))
file.copy(paste(dir_cph_titerAdj,'surv/final/CvR_test.pdf',sep=""),paste(dir_fig_3,'fig_f3e.pdf',sep=""))
file.copy(paste(dir_cph_titerAdj,'surv/final/CvR_compare_test_2.pdf',sep=""),paste(dir_fig_3,'fig_f3f.pdf',sep=""))

# robustness plot
robust_bioph = read.csv(paste(dir_cph_biophysical,'perm/robust.csv',sep=""),row.names = 1)
robust_titer = read.csv(paste(dir_cph_titer,'perm/robust.csv',sep=""),row.names = 1)
robust_titerAdj = read.csv(paste(dir_cph_titerAdj,'perm/robust.csv',sep=""),row.names = 1)

robust_test = as.data.frame(cbind(robust_bioph[1:100,c(2)],robust_titer[1:100,c(2)],robust_titerAdj[1:100,c(2)]))
colnames(robust_test) = c('Fc Array','Titer','Fc Array Adjusted')

robust_perm = as.data.frame(cbind(robust_bioph[101:200,c(2)],robust_titer[101:200,c(2)],robust_titerAdj[101:200,c(2)]))
colnames(robust_perm) = c('Fc Array','Titer','Fc Array Adjusted')

df_test = as.data.frame(cbind(c(rep(1,100),rep(2,100),rep(3,100)),as.vector(as.matrix(robust_test))))
colnames(df_test) = c('label','C_index')
df_test$label = as.factor(df_test$label)

df_all = cbind(robust_test[,1],robust_perm[,1],rep(NA,100),robust_test[,2],robust_perm[,2],rep(NA,100),robust_test[,3],robust_perm[,3])
df_all = as.data.frame(cbind(as.vector(t(1:8 %*% t(rep(1,100)))),as.vector(df_all)))
colnames(df_all) = c('label','C_index')
df_all$label = as.factor(df_all$label)

#---
cdf_perm = ecdf(robust_test[,3])
actual_avg = mean(robust_test[,1])
pval_actual_d1 = 1-cdf_perm(actual_avg)

cdf_perm = ecdf(robust_test[,2])
actual_avg = mean(robust_test[,3])
pval_actual_d2 = 1-cdf_perm(actual_avg)

cdf_perm = ecdf(robust_test[,2])
actual_avg = mean(robust_test[,1])
pval_actual_d3 = 1-cdf_perm(actual_avg)

#---
cdf_perm = ecdf(robust_perm[,1])
actual_avg = mean(robust_test[,1])
pval_perm_d1 = 1-cdf_perm(actual_avg)

cdf_perm = ecdf(robust_perm[,2])
actual_avg = mean(robust_test[,2])
pval_perm_d2 = 1-cdf_perm(actual_avg)

cdf_perm = ecdf(robust_perm[,3])
actual_avg = mean(robust_test[,3])
pval_perm_d3 = 1-cdf_perm(actual_avg)

#---
eff_test_1 = cliff.delta(robust_test[,1],robust_test[,3])
eff_test_2 = cliff.delta(robust_test[,2],robust_test[,3])
eff_test_3 = cliff.delta(robust_test[,1],robust_test[,2])

eff_interp_1 = as.character(eff_test_1$magnitude)
eff_interp_2 = as.character(eff_test_2$magnitude)
eff_interp_3 = as.character(eff_test_3$magnitude)

#---
eff_perm_1 = cliff.delta(robust_test[,1],robust_perm[,1])
eff_perm_2 = cliff.delta(robust_test[,2],robust_perm[,2])
eff_perm_3 = cliff.delta(robust_test[,3],robust_perm[,3])

eff_interp_perm_1 = as.character(eff_perm_1$magnitude)
eff_interp_perm_2 = as.character(eff_perm_2$magnitude)
eff_interp_perm_3 = as.character(eff_perm_3$magnitude)

pdf(paste(dir_fig_3,'fig_f3g.pdf',sep=""))
par(mar=c(12,7,3,1.5))
p = ggplot(df_all,aes(x=label,y=C_index)) + geom_violin(size=1,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Actual','Permuted','','Actual','Permuted','','Actual','Permuted')) + scale_y_continuous(limits = c(0.25,1), breaks=seq(0.25,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=25,colour='black') ,axis.text.x = element_text(angle=25,size=20,colour='black',hjust=1), axis.text.y = element_text(size=17,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='none') + scale_fill_manual(values=c('darkorange3',"#CD660066",'yellow3',"#CDCD00B3",'forestgreen',"#228B22B3")) + xlab("") + ylab('Concordance Index')

p = p + geom_hline(yintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(robust_test[,1],na.rm=T),colour='darkorange3',size=1.5,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(robust_test[,3],na.rm=T),colour='forestgreen',size=1.2,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(robust_test[,2],na.rm=T),colour='yellow3',size=1.2,linetype='dashed',alpha=0.7)

p = p + annotate("segment",x=1,xend=7,y=0.98,yend=0.98,size=2,colour=effect_colors[eff_interp_1])
p = p + annotate("point",shape=24,y=0.98,x=4,size=11,color="black",fill=effect_colors[eff_interp_1])
p = p + annotate("text",x=4,y=0.94,size=6,label=paste('P : ',format(pval_actual_d1,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=4,xend=7,y=0.90,yend=0.90,size=2,colour=effect_colors[eff_interp_2])
p = p + annotate("point",shape=24,y=0.9,x=5.5,size=11,color="black",fill=effect_colors[eff_interp_2])
p = p + annotate("text",x=5.5,y=0.86,size=6,label=paste('P : ',format(pval_actual_d2,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=1,xend=4,y=0.86,yend=0.86,size=2,colour=effect_colors[eff_interp_3])
p = p + annotate("point",shape=24,y=0.86,x=2.5,size=11,color="black",fill=effect_colors[eff_interp_3])
p = p + annotate("text",x=2.5,y=0.82,size=6,label=paste('P : ',format(pval_actual_d3,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=1,xend=2,y=0.35,yend=0.35,size=2,colour=effect_colors[eff_interp_perm_1])
p = p + annotate("point",shape=24,y=0.35,x=1.5,size=9,color="black",fill=effect_colors[eff_interp_perm_1])
p = p + annotate("text",x=1.5,y=0.32,size=6,label=paste('P : ',format(pval_perm_d1,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=4,xend=5,y=0.35,yend=0.35,size=2,colour=effect_colors[eff_interp_perm_2])
p = p + annotate("point",shape=24,y=0.35,x=4.5,size=9,color="black",fill=effect_colors[eff_interp_perm_2])
p = p + annotate("text",x=4.5,y=0.32,size=6,label=paste('P : ',format(pval_perm_d2,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=7,xend=8,y=0.35,yend=0.35,size=2,colour=effect_colors[eff_interp_perm_3])
p = p + annotate("point",shape=24,y=0.35,x=7.5,size=9,color="black",fill=effect_colors[eff_interp_perm_3])
p = p + annotate("text",x=7.5,y=0.32,size=6,label=paste('P : ',format(pval_perm_d3,digits=3,scientific=T),sep=""))
print(p)
dev.off()


# -------------------------------------------
# Sec 05: Figure 4
# -------------------------------------------

dir_fig_4 = paste(dir_res,'fig_4/',sep="")
dir.create(dir_fig_4)

file.copy(paste(dir_class_biophysical,'class/box_pred.pdf',sep=""),paste(dir_fig_4,'fig_f4a.pdf',sep=""))
file.copy(paste(dir_class_biophysical,'class/coeffs_min.pdf',sep=""),paste(dir_fig_4,'fig_f4b.pdf',sep=""))
file.copy(paste(dir_class_biophysical,'class/biplot.pdf',sep=""),paste(dir_fig_4,'fig_f4c.pdf',sep=""))

# robustness plot
robust_bioph = read.csv(paste(dir_class_biophysical,'class/robust.csv',sep=""),row.names = 1)
robust_titer = read.csv(paste(dir_class_titer,'class/robust.csv',sep=""),row.names = 1)
robust_titerAdj = read.csv(paste(dir_class_titerAdj,'class/robust.csv',sep=""),row.names = 1)

robust_test = as.data.frame(cbind(robust_bioph[1:100,c(2)],robust_titer[1:100,c(2)],robust_titerAdj[1:100,c(2)]))
colnames(robust_test) = c('Fc Array','Titer','Fc Array Adjusted')

robust_perm = as.data.frame(cbind(robust_bioph[101:200,c(2)],robust_titer[101:200,c(2)],robust_titerAdj[101:200,c(2)]))
colnames(robust_perm) = c('Fc Array','Titer','Fc Array Adjusted')

df_test = as.data.frame(cbind(c(rep(1,100),rep(2,100),rep(3,100)),as.vector(as.matrix(robust_test))))
colnames(df_test) = c('label','C_index')
df_test$label = as.factor(df_test$label)

df_all = cbind(robust_test[,1],robust_perm[,1],rep(NA,100),robust_test[,2],robust_perm[,2],rep(NA,100),robust_test[,3],robust_perm[,3])
df_all = as.data.frame(cbind(as.vector(t(1:8 %*% t(rep(1,100)))),as.vector(df_all)))
colnames(df_all) = c('label','C_index')
df_all$label = as.factor(df_all$label)

df_all = cbind(robust_test[,1],robust_perm[,1],rep(NA,100),robust_test[,2],robust_perm[,2],rep(NA,100),robust_test[,3],robust_perm[,3])
df_all = as.data.frame(cbind(as.vector(t(1:8 %*% t(rep(1,100)))),as.vector(df_all)))
colnames(df_all) = c('label','C_index')
df_all$label = as.factor(df_all$label)

#---
cdf_perm = ecdf(robust_test[,3])
actual_avg = mean(robust_test[,1])
pval_actual_d1 = 1-cdf_perm(actual_avg)

cdf_perm = ecdf(robust_test[,2])
actual_avg = mean(robust_test[,3])
pval_actual_d2 = 1-cdf_perm(actual_avg)

cdf_perm = ecdf(robust_test[,2])
actual_avg = mean(robust_test[,1])
pval_actual_d3 = 1-cdf_perm(actual_avg)

#---
cdf_perm = ecdf(robust_perm[,1])
actual_avg = mean(robust_test[,1])
pval_perm_d1 = 1-cdf_perm(actual_avg)

cdf_perm = ecdf(robust_perm[,2])
actual_avg = mean(robust_test[,2])
pval_perm_d2 = 1-cdf_perm(actual_avg)

cdf_perm = ecdf(robust_perm[,3])
actual_avg = mean(robust_test[,3])
pval_perm_d3 = 1-cdf_perm(actual_avg)

#---
eff_test_1 = cliff.delta(robust_test[,1],robust_test[,3])
eff_test_2 = cliff.delta(robust_test[,2],robust_test[,3])
eff_test_3 = cliff.delta(robust_test[,1],robust_test[,2])

eff_interp_1 = as.character(eff_test_1$magnitude)
eff_interp_2 = as.character(eff_test_2$magnitude)
eff_interp_3 = as.character(eff_test_3$magnitude)

#---
eff_perm_1 = cliff.delta(robust_test[,1],robust_perm[,1])
eff_perm_2 = cliff.delta(robust_test[,2],robust_perm[,2])
eff_perm_3 = cliff.delta(robust_test[,3],robust_perm[,3])

eff_interp_perm_1 = as.character(eff_perm_1$magnitude)
eff_interp_perm_2 = as.character(eff_perm_2$magnitude)
eff_interp_perm_3 = as.character(eff_perm_3$magnitude)

pdf(paste(dir_fig_4,'fig_f4d.pdf',sep=""))
par(mar=c(12,7,3,1.5))
p = ggplot(df_all,aes(x=label,y=C_index)) + geom_violin(size=1,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Actual','Permuted','','Actual','Permuted','','Actual','Permuted')) + scale_y_continuous(limits = c(0.15,1.1), breaks=seq(0.1,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=25,colour='black') ,axis.text.x = element_text(angle=25,size=20,colour='black',hjust=1), axis.text.y = element_text(size=17,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='none') + scale_fill_manual(values=c('darkorange3',"#CD660066",'yellow3',"#CDCD00B3",'forestgreen',"#228B22B3")) + xlab("") + ylab('Balanced Accuracy')

p = p + geom_hline(yintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(robust_test[,1],na.rm=T),colour='darkorange3',size=1.5,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(robust_test[,3],na.rm=T),colour='forestgreen',size=1.2,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(robust_test[,2],na.rm=T),colour='yellow3',size=1.2,linetype='dashed',alpha=0.7)

p = p + annotate("segment",x=1,xend=7,y=1.08,yend=1.08,size=2,colour=effect_colors[eff_interp_1])
p = p + annotate("point",shape=24,y=1.08,x=4,size=11,color="black",fill=effect_colors[eff_interp_1])
p = p + annotate("text",x=4,y=1.04,size=6,label=paste('P : ',format(pval_actual_d1,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=4,xend=7,y=1.0,yend=1.0,size=2,colour=effect_colors[eff_interp_2])
p = p + annotate("point",shape=24,y=1.0,x=5.5,size=11,color="black",fill=effect_colors[eff_interp_2])
p = p + annotate("text",x=5.5,y=0.96,size=6,label=paste('P : ',format(pval_actual_d2,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=1,xend=4,y=0.98,yend=0.98,size=2,colour=effect_colors[eff_interp_3])
p = p + annotate("point",shape=24,y=0.98,x=2.5,size=11,color="black",fill=effect_colors[eff_interp_3])
p = p + annotate("text",x=2.5,y=0.94,size=6,label=paste('P : ',format(pval_actual_d3,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=1,xend=2,y=0.27,yend=0.27,size=2,colour=effect_colors[eff_interp_perm_1])
p = p + annotate("point",shape=24,y=0.27,x=1.5,size=9,color="black",fill=effect_colors[eff_interp_perm_1])
p = p + annotate("text",x=1.5,y=0.23,size=6,label=paste('P : ',format(pval_perm_d1,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=4,xend=5,y=0.27,yend=0.27,size=2,colour=effect_colors[eff_interp_perm_2])
p = p + annotate("point",shape=24,y=0.27,x=4.5,size=9,color="black",fill=effect_colors[eff_interp_perm_2])
p = p + annotate("text",x=4.5,y=0.23,size=6,label=paste('P : ',format(pval_perm_d2,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=7,xend=8,y=0.27,yend=0.27,size=2,colour=effect_colors[eff_interp_perm_3])
p = p + annotate("point",shape=24,y=0.27,x=7.5,size=9,color="black",fill=effect_colors[eff_interp_perm_3])
p = p + annotate("text",x=7.5,y=0.23,size=6,label=paste('P : ',format(pval_perm_d3,digits=3,scientific=T),sep=""))
print(p)
dev.off()

# -------------------------------------------
# Sec 06: Figure S1
# -------------------------------------------

# -------------------------------------------
# Sec 07: Figure S2
# -------------------------------------------

# -------------------------------------------
# Sec 08: Figure S3
# -------------------------------------------

dir_fig_s3 = paste(dir_res,'fig_s3/',sep="")
dir.create(dir_fig_s3)

file.copy(paste(dir_tcell_relation,'IFNg_relation.pdf',sep=""),paste(dir_fig_s3,'fig_s3a.pdf',sep=""))
file.copy(paste(dir_tcell_relation,'ADCC_relation.pdf',sep=""),paste(dir_fig_s3,'fig_s3b.pdf',sep=""))
file.copy(paste(dir_tcell_relation,'separation_by_risk.pdf',sep=""),paste(dir_fig_s3,'fig_s3c.pdf',sep=""))

# -------------------------------------------
# Sec 09: Figure S4
# -------------------------------------------

# -------------------------------------------
# Sec 10: Figure S5
# -------------------------------------------

dir_fig_s5 = paste(dir_res,'fig_s5/',sep="")
dir.create(dir_fig_s5)

file.copy(paste(dir_multi_class,'class/confusion.pdf',sep=""),paste(dir_fig_s5,'fig_s5a.pdf',sep=""))
file.copy(paste(dir_multi_class,'class/box_odds_gg.pdf',sep=""),paste(dir_fig_s5,'fig_s5b.pdf',sep=""))
file.copy(paste(dir_multi_class,'class/robust_test_v.pdf',sep=""),paste(dir_fig_s5,'fig_s5c.pdf',sep=""))
file.copy(paste(dir_multi_class,'class/coeffs_min.pdf',sep=""),paste(dir_fig_s5,'fig_s5d.pdf',sep=""))

cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)
