# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectus T2
# Decription  : Build a final model using the input features
# Cite        : TBD
# ******************************************************************************

coxSurvivalFinalModel = function(feats,subj_original,seedMe,scolors,group_colors,km_colors,num_folds,plots,dir_surv){
  
  plot_font = 'Helvetica'
  
  survObject = Surv((subj_original[,'Challenges']+subj_original[,'censor']-1),subj_original[,'censor'])
  
  # -------------
  # Standardize
  # -------------
  # if subjects have NA in features, set those features to 0
  na_idx = which(is.na(feats),arr.ind=TRUE)
  feats_scaled = scale(data.matrix(feats),center=TRUE,scale=TRUE)
  if(length(na_idx) !=0)
    feats_scaled[na_idx] = 0
  
  feats_scaled = data.frame(feats_scaled)
  
  group_ID = unique(subj_original[,'groupID'])
  
  # learn a final cox model using all subjects
  fit_cox_all = coxph(survObject ~ .,data=feats_scaled,control=coxph.control(iter.max=100))
  km_all = survfit(survObject ~ 1)
  
  # predict KM curves for each group using the final cox model
  pdf(paste(dir_surv,'surv_cox_group.pdf',sep=""),family=plot_font)
  plot(survfit(fit_cox_all,newdata=colMeans(feats_scaled[which(subj_original[,'groupID']==1),])),mark.time=FALSE,conf.int=FALSE,col=km_colors[1],lwd=3,xlab='Challenges',ylab='uninfected ratio')
  lines(survfit(fit_cox_all,newdata=colMeans(feats_scaled[which(subj_original[,'groupID']==2),])),mark.time=FALSE,conf.int=FALSE,col=km_colors[2],lwd=3)
  dev.off()
  
  max_chall = max(subj_original[,'Challenges'] + subj_original[,'censor'] - 1)
  
  #-----------------------------------------------
  # Quantify comparison between actaul and predicted KM curves
  #-----------------------------------------------
  
  # -------------
  # actual survival probabilities
  km_actual = survfit(survObject ~ strata(subj_original[,'groupID']))
  event_size = cumsum(km_actual$strata)
  start_idx = c(1,1+event_size)
  
  surv_actual_grp = as.list(rep(NA,length(group_ID)))
  names(surv_actual_grp) = paste('grp_',group_ID,sep="")
  
  for(grpIdx in 1:length(group_ID)){
    
    time_grp = km_actual$time[start_idx[[grpIdx]]:event_size[[grpIdx]]]
    surv_grp = km_actual$surv[start_idx[[grpIdx]]:event_size[[grpIdx]]]
    
    surv_actual_grp[[grpIdx]] = extractProbabilityFromKM(surv_grp,time_grp,max_chall)
    
  }
  
  # -------------
  # predicted survival probabilities
  surv_pred_grp = as.list(rep(NA,length(group_ID)))
  names(surv_pred_grp) = paste('grp_',group_ID,sep="")
  
  for(grpIdx in 1:length(group_ID)){
    
    km_pred_grp = survfit(fit_cox_all,newdata=colMeans(feats_scaled[which(subj_original[,'groupID']==group_ID[grpIdx]),]))
    surv_pred_grp[[grpIdx]] = extractProbabilityFromKM(km_pred_grp$surv,km_pred_grp$time,max_chall)
    
  }
  
  # -------------
  # plot comparison between actual and predicted KM curves
  
  events_pred_grp = as.list(rep(NA,length(group_ID)))
  names(events_pred_grp) = paste('grp_',group_ID,sep="")
  
  table_pred_grp = as.list(rep(NA,length(group_ID)))
  names(table_pred_grp) = paste('grp_',group_ID,sep="")
  
  surv_grp = as.list(rep(NA,length(group_ID)))
  names(surv_grp) = paste('grp_',group_ID,sep="")
  
  pval_grp = as.list(rep(NA,length(group_ID)))
  names(pval_grp) = paste('grp_',group_ID,sep="")
  
  for(grpIdx in 1:length(group_ID)){
    
    grp_ID = which(subj_original[,'groupID']==group_ID[grpIdx])
    num_samples = length(grp_ID)
    flag = c(((1:num_samples)*0)+1,((1:num_samples)*0))
    
    events_pred_grp[[grpIdx]] = convertKMtoEvents(surv_pred_grp[[grpIdx]],max_chall,num_samples)
    write.csv(events_pred_grp[[grpIdx]]$event_table,file=paste(dir_surv,names(events_pred_grp)[grpIdx],'_pred.csv',sep=""),row.names = T)
    
    table_pred_grp[[grpIdx]] = rbind(cbind(subj_original[grp_ID,'Challenges'],subj_original[grp_ID,'censor']),events_pred_grp[[grpIdx]]$event_table)
    surv_grp[[grpIdx]] = survdiff(Surv(table_pred_grp[[grpIdx]][,'Challenges'],table_pred_grp[[grpIdx]][,'censor']) ~ flag)
    pval_grp[[grpIdx]] = round(1 - pchisq(surv_grp[[grpIdx]]$chisq,1),2)
    
  }
  
  #-----------------------------------------------
  # Write output
  
  label = NULL
  surv_actual_all = NULL
  surv_pred_all = NULL
  for(grpIdx in 1:length(group_ID)){
    
    label = c(label,rep(grpIdx,max_chall+1))
    surv_actual_all = c(surv_actual_all,c(1,surv_actual_grp[[grpIdx]]))
    surv_pred_all = c(surv_pred_all,c(1,surv_pred_grp[[grpIdx]]))
    
  }
  
  km_points = cbind(label,cbind(surv_actual_all,surv_pred_all))
  pvals = unlist(pval_grp)
  
  colnames(km_points) = c('label','Actual','Pred')
  
  write.csv(rbind(km_points,pvals),file=paste(dir_surv,'km_points.csv',sep=""),row.names = F)
  
  len_table = length(surv_pred_all)
  
  Challenges = rep(c(0,seq(1:max_chall)),length(group_ID))
  km_points = cbind(Challenges,km_points)
  new = rbind(as.matrix(km_points[,c(1,2,3)]),as.matrix(km_points[,c(1,2,4)]))
  type = c(rep(1,len_table),rep(2,len_table))
  km_points = as.data.frame(cbind(type,new))
  km_points$type = as.factor(km_points$type)
  km_points$label = as.factor(km_points$label)
  
  pdf(paste(dir_surv,'km_compare.pdf',sep=""),family=plot_font)
  par(mar=c(12,4,2,0.5))
  p = ggplot(km_points,aes(x=Challenges, y=Actual, color=label)) + geom_step(size=1.2,aes(linetype=type)) + scale_x_discrete(breaks=0:max_chall,limits=0:max_chall,labels=c(as.character(0:max_chall)),expand = c(0,0.2)) + scale_y_continuous(breaks=seq(0,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=25,colour='black'), axis.title.y = element_text(size=25,colour='black') ,axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=20,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='bottom') + scale_colour_manual(values=km_colors) + scale_linetype_manual(values=c('dashed','solid')) + xlab('Challenges') + ylab('Survival Probability')
  #p = p + annotate("text",x = c(8,3), y = c(0.21,0.96),label=paste('p: ',c(pvals[[1]],pvals[[2]]),sep=""),color=km_colors,size=7)
  p = p + expand_limits(x = 0, y = 0)
  print(p)
  dev.off()
  
  #-----------------------------------------------
  # plot predictor coefficients
  
  coeffs = fit_cox_all$coefficients
  coeffs = exp(coeffs[which(coeffs!=0)])
  
  #
  # statistical summaries
  cox_summary = summary(fit_cox_all)
  cox_signif = cox_summary$coefficients[,5]
  cox_signif_sym = character(length(cox_signif))
  names(cox_signif_sym) = names(coeffs)
  for(featID in 1:length(cox_signif)){
    
    if(cox_signif[featID]<0.001){
      cox_signif_sym[featID] = '***'
    }else if(cox_signif[featID]<0.01){
      cox_signif_sym[featID] = '**'
    }else if(cox_signif[featID]<0.05){
      cox_signif_sym[featID] = '*'
    }else if(cox_signif[featID]<0.1){
      cox_signif_sym[featID] = '-'
    }else{
      cox_signif_sym[featID] = '-'
    }
    
  }
  
  # -------------
  # make heatmap of final feature set
  
  coeff_order = order(coeffs,decreasing=TRUE)
  feat_data = data.matrix(feats_scaled[,coeff_order])
  feat_data[feat_data>3] = 3; feat_data[feat_data < -3] = -3
  fr = ceiling(max(abs(min(feat_data,na.rm=TRUE)), max(feat_data,na.rm=TRUE)))
  fbreaks = seq(-fr,fr,0.1)
  
  coeff_sign = matrix(nrow=ncol(feats_scaled), ncol=1, dimnames=list(colnames(feats_scaled),c('sign')))
  coeff_sign[which(coeffs[coeff_order]<=1)] = -0.1
  coeff_sign[which(coeffs[coeff_order]>1)] = 0.1
  
  wcolors = matrix(nrow=ncol(feat_data), ncol=1, dimnames=list(colnames(feat_data),c('weights')))
  wcolors[which(coeffs[coeff_order]<=1),'weights'] = 'darkgoldenrod2'
  wcolors[which(coeffs[coeff_order]>1),'weights'] = 'khaki4'
  
  coeffs_all = coeffs[coeff_order]
  pdf(paste(dir_surv,'surv_cox_all.pdf',sep=""),family=plot_font)
  if(length(coeffs_all != 0)){
    par(mar=c(12,7,3,1.5))
    yrange = range(log(coeffs_all))
    bp = barplot(log(coeffs_all),main=c(""),las=2,names.arg=colnames(feat_data),col=wcolors,cex.names=1,width=0.75,xlim=c(0,length(coeffs_all)+1),ylim=c(-1,1),ylab="",cex.axis=2)
    mtext(expression(paste('Coefficient (',beta,')',sep="")), side=2, line=4, cex=3)
    text(x=bp,y=log(coeffs_all)+coeff_sign,labels=cox_signif_sym[coeff_order],cex=4.5,xpd=TRUE)
    abline(h=0,lty=2,lwd=2)
  }else{
    plot(c(0,1),c(0,1),type='n',xaxt='n',yaxt='n',xlab='',ylab=''); text(0.5,0.5,'(empty model)')
  }
  dev.off()
  
  pdf(paste(dir_surv,'final_feat_selection_sorted.pdf',sep=""),family=plot_font)
  chall.sort <- sort(subj_original[,'Challenges'], index.return=TRUE)
  heatmap.4(feat_data[chall.sort$ix,], col=bluered, scale='none', trace='none', cexRow=0.5, cexCol=1.2, margin=c(14,5), breaks=fbreaks, symkey=FALSE,dendrogram='none',RowSideColors=scolors[chall.sort$ix,],NumRowSideColors=4,Rowv=FALSE,ColSideColors=wcolors, NumColSideColors=1,Colv=FALSE,na.color='black',lmat=rbind(c(6,0,5,0),c(0,0,2,0),c(4,1,3,0)), lhei=c(1.5,0.3,7.0),lwid=c(0.2,0.06,0.2,0.3))
  dev.off()
  
  pdf(paste(dir_surv,'final_feat_selection_raw.pdf',sep=""))
  heatmap.4(feat_data, col=bluered, scale='none', trace='none', cexRow=0.5, cexCol=1, margin=c(14,5), breaks=fbreaks, symkey=FALSE, dendrogram='none',RowSideColors=scolors,NumRowSideColors=4,Rowv=FALSE,ColSideColors=wcolors, NumColSideColors=1,Colv=FALSE,na.color='black',lmat=rbind(c(6,0,5,0),c(0,0,2,0),c(4,1,3,0)), lhei=c(2,0.2,6.0),lwid=c(0.2,0.06,0.2,0.3))
  dev.off()
  
  #---------------------------------------------
  # Cross Validation
  
  dir_cv = paste(dir_surv,"cv/",sep="")
  dir.create(dir_cv)
  
  cat(num_folds,'Fold CV\n')
  
  val_subjects = subj_original
  num_samples = nrow(val_subjects)
  
  # -------------
  # collect predictions
  cox_pred_lp = array(0,dim=c(num_samples,1))
  cox_pred_lp_train = matrix(NA,nrow=num_samples,ncol=num_folds)
  
  set.seed(seedMe)
  
  folds_list = createFolds(val_subjects[,'groupID'],num_folds,list=FALSE)
  
  overall_mean = colMeans(feats,na.rm=TRUE)
  
  for(foldIdx in 1:num_folds){
    
    dir_fold = paste(dir_cv,"fold_",foldIdx,"/",sep="")
    dir.create(dir_fold)
    
    train_idx = which(folds_list != foldIdx)
    test_idx = which(folds_list == foldIdx)
    
    # -------------
    # prepare trainng data
    cv_survObject = Surv((val_subjects[train_idx,'Challenges']+val_subjects[train_idx,'censor']-1),val_subjects[train_idx,'censor'])
    
    cv_center = preProcess(feats[train_idx,],method=c('center','scale'),na.rm=TRUE)
    cv_feats = predict(cv_center,feats[train_idx,])
    cv_test_feats = predict(cv_center,feats[test_idx,])
    
    overall_mean_sc = predict(cv_center,t(replicate(length(test_idx),overall_mean)))
    
    cv_feats[which(is.na(cv_feats),arr.ind=TRUE)] = 0
    cv_test_feats[which(is.na(cv_test_feats),arr.ind=TRUE)] = 0
    
    # -------------
    # learn a cox model
    fit_cox_cv = coxph(cv_survObject ~ .,data=cv_feats,control=coxph.control(iter.max=100))
    
    # -------------
    # predict on train data
    cox_pred_train_cv = predict(fit_cox_cv,newdata=(cv_feats-predict(cv_center,t(replicate(length(train_idx),overall_mean)))),type="lp")
    cox_pred_lp_train[train_idx,foldIdx] = cox_pred_train_cv
    
    # -------------
    # predict on test data
    cox_pred_cv = predict(fit_cox_cv,newdata=(cv_test_feats-overall_mean_sc),type="lp")
    cox_pred_lp[test_idx,1] = cox_pred_cv
    
  }
  rm(.Random.seed,envir=globalenv())
  
  #cox_pred_lp[cox_pred_lp>1.9] = 1.9;
  
  # Observed vs Predicted
  OvP = as.data.frame(cbind(val_subjects[,'Challenges'],cox_pred_lp,rowSums(cox_pred_lp_train,na.rm=TRUE)/(num_folds-1),val_subjects[,'groupID'],val_subjects[,'group_all']))
  colnames(OvP) = c("Challenges","Risk","Risk_Train","Group","Group_all")
  
  # -------------
  # test set
  cat(rep("#",10)," Test Set ",rep("#",11),"\n",sep="")
  ppc_test = round(polyserial(OvP[,"Risk"],OvP[,"Challenges"]),2)
  cindex_test = concordance.index(OvP$Risk,(val_subjects[,'Challenges']+val_subjects[,'censor']-1),val_subjects[,'censor'],method="noether")
  cat('#Overall Polyserial : ',ppc_test,'\n')
  cindex_ov = paste(round(cindex_test$c.index,2),' p: ',format(cindex_test$p.value,scientific=T))
  cat('#Overall C-index is : ',cindex_ov,'\n')
  
  # -------------
  # train set
  cat(rep("#",10)," Train Set ",rep("#",10),"\n",sep="")
  ppc_train = round(polyserial(OvP[,"Risk_Train"],OvP[,"Challenges"]),2)
  cindex_train = concordance.index(OvP$Risk_Train,(val_subjects[,'Challenges']+val_subjects[,'censor']-1),val_subjects[,'censor'],method="noether")
  cat('#Overall Polyserial : ',ppc_train,'\n')
  cat('Overall C-index is : ',round(cindex_train$c.index,2),' p : ',cindex_train$p.value,'\n')
  
  model = lm(Risk_Train~Risk,OvP)
  pcc_val = round(sqrt(summary(model)$r.squared),2)
  cat("Train vs Test pcc = ",pcc_val,"\n")
  
  write.csv(OvP,file=paste(dir_surv,'ovp.csv',sep=""),row.names = T,col.names = T)
  
  if(plots){
    
    pdf(paste(dir_surv,"CvR_test.pdf",sep=""),family=plot_font)
    OvP$Group = as.factor(OvP$Group)
    OvP$Group_all = as.factor(OvP$Group_all)
    
    p = ggplot(OvP, aes(x=Challenges, y=Risk, color=Group_all)) + geom_smooth(method="lm",se=F,formula=y~log(x),inherit.aes=F,aes(x=Challenges,y=Risk),colour='black',size=1.5,linetype='dashed') + geom_point(aes(size=Group)) + scale_x_discrete(breaks=1:(max_chall+1),limits=1:(max_chall+1),labels=c(as.character(1:max_chall),'UI')) + ylab('Predicted Risk') + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=25,colour='black'), axis.title.y = element_text(size=25,colour='black') ,axis.text.x = element_text(size=15,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_line(colour='gray65',size=0.3,linetype = 'dashed'), panel.grid.minor.y = element_blank(), aspect.ratio=1, panel.grid.major.y = element_blank(), legend.position='none') + scale_colour_manual(name='Group\nC-index',values=group_colors) + xlab('Challenges') + scale_size_manual(name='Overall\nC-index',values=rep(3,length(group_ID)),labels=c(cindex_ov))
    p = p + geom_hline(yintercept=0,colour='black',size=0.78,linetype='dashed',alpha=0.7)
    print(p)
    dev.off()
    
    # Split CvR group-wise
    pdf(paste(dir_surv,"CvR_compare_test_2.pdf",sep=""),family=plot_font)
    yrange = range(OvP$Risk,na.rm=T)
    
    
    grp_1 = which(OvP[,'Group']==1)
    grp_2 = which(OvP[,'Group']==2)
    
    diff_g1g2 = wilcox.test(OvP[grp_1,'Risk'],OvP[grp_2,'Risk'],alternative="two.sided")
    
    p = ggplot(OvP, aes(x=Group, y=Risk, color=Group_all)) + geom_boxplot(width=0.4,notch = F,outlier.shape = NA, na.rm=T, size=1,colour="black") + geom_point(position = position_jitter(w=0.1),size=4) + scale_x_discrete(labels=names(group_id)) + scale_y_continuous(limits=c(yrange[1],yrange[2]*1.7)) + ylab('Predicted Risk') + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=25,colour='black'), axis.title.y = element_text(size=25,colour='black') ,axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),aspect.ratio=1, legend.position='none') + xlab('') + scale_colour_manual(values=group_colors)
    p = p + geom_hline(yintercept=0,colour='black',size=0.78,linetype='dashed',alpha=0.7)
    
    p = p + annotate("segment",x=1,xend=2,y=yrange[2]*1.1,yend=yrange[2]*1.1,size=1.5,colour='black')
    p = p + annotate("text",x=1.5,y=yrange[2]*1.18,size=5,label=paste('P : ',format(diff_g1g2$p.value,digits=3,scientific=T)))
    
    print(p)
    dev.off()
    
    pdf(paste(dir_surv,"CvR_train.pdf",sep=""))
    p = ggplot(OvP, aes(x=Challenges, y=Risk_Train, color=Group_all)) + geom_point(size=3) + ggtitle(paste("Challenges vs Predicted Risk (polyserial : ",ppc_train,")",sep="")) + scale_colour_manual(name='Groupwise Coeff\nPolyserial',values=group_colors)
    print(p)
    dev.off()
    
    pdf(paste(dir_surv,"RvR.pdf",sep=""))
    p = ggplot(OvP, aes(x=Risk_Train, y=Risk)) + geom_point() + ggtitle(paste("Predicted Risk - Train vs Test (Regr coeff : ",pcc_val,")",sep="")) + coord_fixed(ratio=1) + geom_smooth(method="lm",fullrange=TRUE,se=FALSE)
    print(p)
    dev.off()
    
  }
  
  return(list('ppc_train'=ppc_train,'cindex_train'=round(cindex_train$c.index,2),'ppc_test'=ppc_test,'cindex_test'=round(cindex_test$c.index,2),'pcc_val'=pcc_val,'coeffs'=coeffs_all))
  
}