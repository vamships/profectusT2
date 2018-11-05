# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectus T2
# Decription  : Performs survival analysis using the input features and
#               backward feature selection
# Cite        : TBD
# ******************************************************************************

coxSurvivalWithBackSearch = function(feats,subj_original,selected_feat_idx,group_colors,num_folds,stop_limit,plots,dir_surv){
  
  # -------------
  # Standardize
  # -------------
  # if subjects have NA in features, set those features to 0
  na_idx = which(is.na(feats),arr.ind=TRUE)
  feats_scaled = scale(data.matrix(feats),center=TRUE,scale=TRUE)
  if(length(na_idx) !=0)
    feats_scaled[na_idx] = 0
  
  feats_scaled = data.frame(feats_scaled)
  
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
  
  # -------------
  #stratified sampling based on groupID
  folds_list = createFolds(val_subjects[,'groupID'],num_folds,list=FALSE)
  
  fold_feats_count = matrix(0,nrow=num_folds,ncol=length(selected_feat_idx))
  colnames(fold_feats_count) = names(selected_feat_idx)
  
  # -------------
  # mean feature vector over all samples
  overall_mean = colMeans(feats,na.rm=TRUE)
  
  for(foldIdx in 1:num_folds){
    
    dir_fold = paste(dir_cv,"fold_",foldIdx,"/",sep="")
    
    if(plots){
      
      dir.create(dir_fold)
      
    }
    
    
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
    #select features for the fold
    fold_feat_idx = doCoxBackSearch(cv_feats,cv_survObject,selected_feat_idx,stop_limit,plots,dir_fold)
    sub_select_idx = which(selected_feat_idx %in% fold_feat_idx)
    
    fold_feats_count[foldIdx,sub_select_idx] = 1
    
    # -------------
    # prepare test data
    cv_feats = cv_feats[,sub_select_idx]
    cv_test_feats = cv_test_feats[,sub_select_idx]
    overall_mean_sc = overall_mean_sc[,sub_select_idx]
    
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
  
  # Observed vs Predicted
  OvP = as.data.frame(cbind(val_subjects[,'Challenges'],cox_pred_lp,rowSums(cox_pred_lp_train,na.rm=TRUE)/(num_folds-1),val_subjects[,'groupID']))
  colnames(OvP) = c("Challenges","Risk","Risk_Train","Group")
  
  # -------------
  # Test set performance
  cat(rep("#",10)," Test Set ",rep("#",11),"\n",sep="")
  ppc_test = round(polyserial(OvP[,"Risk"],OvP[,"Challenges"]),2)
  cindex_test = concordance.index(OvP$Risk,(val_subjects[,'Challenges']+val_subjects[,'censor']-1),val_subjects[,'censor'],method="noether")
  cat('#Overall C-index is : ',round(cindex_test$c.index,2),' p : ',cindex_test$p.value,'\n')
  cat('#Overall Polyserial : ',ppc_test,'\n')
  
  # -------------
  # Training set performance
  cat(rep("#",10)," Train Set ",rep("#",10),"\n",sep="")
  ppc_train = round(polyserial(OvP[,"Risk_Train"],OvP[,"Challenges"]),2)
  cindex_train = concordance.index(OvP$Risk_Train,(val_subjects[,'Challenges']+val_subjects[,'censor']-1),val_subjects[,'censor'],method="noether")
  cat('#Overall C-index is : ',round(cindex_train$c.index,2),' p : ',cindex_train$p.value,'\n')
  cat('#Overall Polyserial : ',ppc_train,'\n')
  
  # -------------
  # correlation between test and train set predictions
  model = lm(Risk_Train~Risk,OvP)
  pcc_val = round(sqrt(summary(model)$r.squared),2)
  cat("Train vs Test pcc = ",pcc_val,"\n")
  
  if(plots){
    
    pdf(paste(dir_surv,"CvR_test.pdf",sep=""))
    OvP$Group = as.factor(OvP$Group)
    p = ggplot(OvP, aes(x=Challenges, y=Risk, color=Group)) + geom_point(size=3) + ggtitle(paste("Challenges vs Predicted Risk (C-index : ",cindex_test,")",sep="")) + scale_colour_manual(name='Groupwise Coeff\nPolyserial',values=group_colors)
    print(p)
    dev.off()
    
    pdf(paste(dir_surv,"CvR_train.pdf",sep=""))
    p = ggplot(OvP, aes(x=Challenges, y=Risk_Train, color=Group)) + geom_point(size=3) + ggtitle(paste("Challenges vs Predicted Risk (C-index : ",cindex_train,")",sep="")) + scale_colour_manual(name='Groupwise Coeff\nPolyserial',values=group_colors)
    print(p)
    dev.off()
    
    pdf(paste(dir_surv,"RvR.pdf",sep=""))
    p = ggplot(OvP, aes(x=Risk_Train, y=Risk)) + geom_point() + ggtitle(paste("Predicted Risk - Train vs Test (Regr coeff : ",pcc_val,")",sep="")) + coord_fixed(ratio=1) + geom_smooth(method="lm",fullrange=TRUE,se=FALSE)
    print(p)
    dev.off()
    
  }
  
  return(list('feat_counts'=fold_feats_count,'poly_test'=ppc_test,'cindex_test'=round(cindex_test$c.index,2),'poly_train'=ppc_train,'cindex_train'=round(cindex_train$c.index,2),'regr_coeff'=pcc_val))
  
}