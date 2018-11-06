# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectus T2
# Description : Find relation between the risk predictions and cellular functions previously published
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

dir_res = paste('results_func_relation/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

# -------------------------------------------
# Sec 02: Data
# -------------------------------------------

luminex = read.csv('in/luminex_filtered.csv', header=TRUE, row.names=1)

subjects = read.csv('in/subjects.csv', header=TRUE, row.names=1)

cell_response = read.csv('in/func_response.csv', header=TRUE, row.names=1)

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
# Sec 04: Compare biophysical measurement based risk prediction to functional response
# -------------------------------------------

dir_final = paste("results_coxPred_biophysical/surv/final/",sep="")

sel_IFNg = cell_response[,1]
q3 = summary(sel_IFNg)['3rd Qu.']

sel_adcc = cell_response[,2]

pred_final = read.csv(paste(dir_final,'ovp.csv',sep=""),row.names = 1,header=T)
risk_test = pred_final[,'Risk']

# -------------------------------------------
# Survival separation by risk and IFNg
# -------------------------------------------

risk_pos_idx = which(risk_test > 0)
risk_neg_idx = which(risk_test < 0)

# those that have negative risk and more than 3rd quartile
risk_neg_high_idx = intersect(risk_neg_idx,which(sel_IFNg>=q3))

# those that have negative risk and less than 3rd quartile
g23_label = (sel_IFNg[risk_neg_idx]<q3)*1

# the remaining subjects
risk_neg_low_idx = setdiff(risk_neg_idx,risk_neg_high_idx)
g13_idx = c(risk_pos_idx,risk_neg_low_idx)
g13_label = c(rep(0,length(risk_pos_idx)),rep(1,length(risk_neg_low_idx)))

pdf(paste(dir_res,'separation_by_risk.pdf',sep=""))
surv = Surv(subjects[,'Challenges']+subjects[,'censor']-1, subjects[,'censor'])

surv_g13 = Surv(subjects[g13_idx,'Challenges']+subjects[g13_idx,'censor']-1, subjects[g13_idx,'censor'])
km_diff_g13 = survdiff(surv_g13~g13_label)
pval_g13 = round(1 - pchisq(km_diff_g13$chisq,1),3)

surv_g23 = Surv(subjects[risk_neg_idx,'Challenges']+subjects[risk_neg_idx,'censor']-1, subjects[risk_neg_idx,'censor'])
km_diff_g23 = survdiff(surv_g23~g23_label)
pval_g23 = round(1 - pchisq(km_diff_g23$chisq,1),3)

plot(survfit(surv ~ interaction(risk_test>0, risk_test<0 & sel_IFNg>=q3, risk_test<0 & sel_IFNg<q3)), mark.time=FALSE, col=c('blue','green4','red'), lwd=4, main=paste('Log rank p(blue,red): ',pval_g13,' p(green,red): ',pval_g23,sep=""),bty='l',cex.axis=2)
dev.off()

# -------------------------------------------
# IFNg relation
# -------------------------------------------

risk_ind = rep(1,nrow(subjects))
risk_ind[risk_neg_idx] = 2
df = as.data.frame(cbind(cell_response[,1],subjects[,c(1,2)],risk_ind))
colnames(df) = c('vf','groupID','Challenges','risk')
df$groupID = as.factor(df$groupID)
df$risk = factor(df$risk)

pdf(paste(dir_res,"IFNg_relation.pdf",sep=""))
p1 = ggplot(df,aes(x=Challenges,y=vf),color="black") + geom_point(size=4,aes(shape=risk),fill="black") + ylab('MFI') + xlab('Challenges') + ggtitle(paste("IFNg \nChallenges vs Feat",sep="")) + scale_x_discrete(breaks=1:11,limits=1:11,labels=c(as.character(1:10),'UI')) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=25,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=15,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_line(colour='gray65',size=0.3,linetype = 'dashed'), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='bottom') + scale_shape_manual(values=c(5,23)) + scale_colour_manual(name='Groupwise Coeff\nSpearman',values=group_colors)
p1 = p1 + geom_hline(yintercept=180.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
print(p1)
dev.off()

# -------------------------------------------
# ADCC relation
# -------------------------------------------

risk_ind = rep(1,nrow(subjects))
risk_ind[risk_neg_idx] = 2
df = as.data.frame(cbind(cell_response[,2],subjects[,c(1,2)],risk_ind))
colnames(df) = c('vf','groupID','Challenges','risk')
df["p3.3","vf"] = -50
df["p2.8","vf"] = -50
df["p3.4","vf"] = -100
df["p3.5","vf"] = -150
df$groupID = as.factor(df$groupID)
df$risk = factor(df$risk)

pdf(paste(dir_res,"ADCC_relation.pdf",sep=""))
p1 = ggplot(df,aes(x=Challenges,y=vf),color="black") + geom_point(size=4,aes(shape=risk),fill="black") + ylab('MFI') + xlab('Challenges') + ggtitle(paste("ADCC\nChallenges vs Feat",sep="")) + scale_x_discrete(breaks=1:11,limits=1:11,labels=c(as.character(1:10),'UI')) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=25,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=15,colour='black'), axis.text.y = element_text(size=12,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_line(colour='gray65',size=0.3,linetype = 'dashed'), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='bottom') + scale_shape_manual(values=c(5,23)) + scale_colour_manual(name='Groupwise Coeff\nSpearman',values=group_colors)
p1 = p1 + geom_hline(yintercept=0,colour='black',size=0.78,linetype='dashed',alpha=0.7)
print(p1)
dev.off()


cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)
