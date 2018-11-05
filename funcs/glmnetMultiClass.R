glmnetMultiClass = function(feats,predVec,weights,numFeat,grpType,intc,alphas,cvFolds,repeatRun){
  
  augMatrix = cbind(feats,predVec)
  
  rmList = which(is.na(predVec))
  
  if(length(rmList)!=0){
    
    if(cvFolds==nrow(augMatrix))
      cvFolds = cvFolds - length(rmList)
    
    augMatrix = augMatrix[-rmList,]
  }
  
  cat("removed : ",length(rmList)," subjects.\n")
  cat("num folds : ",cvFolds,"\n")
  
  feats = augMatrix[,1:numFeat]
  func = augMatrix[,(numFeat+1)]
  
  # Cross-validation
  alpha_mse = matrix(numeric(1),length(alphas))
  
  foldID = createFolds(func,k=cvFolds,list=F)
  
  cv_glm = cv.glmnet(feats,func,nfolds=cvFolds,foldid=foldID,family="multinomial",standardize=FALSE,weights=weights,alpha=alphas[1], type.multinomial = grpType,type.measure="class",keep=TRUE, intercept=intc)
  
  alphaID = 1
  alpha_mse[1] = min(cv_glm$cvm)
  best_alpha_mse = min(cv_glm$cvm)
  best_alpha = alphas[1]
  foldID = cv_glm$foldid
  
  if(length(alphas)>1){
    
    for(alphax in alphas[2:length(alphas)]){
      
      alphaID = alphaID + 1
      
      cv_glm = cv.glmnet(feats,func,nfolds=cvFolds,foldid=foldID,family="multinomial",standardize=FALSE,weights=weights,alpha=alphax, type.multinomial = grpType,type.measure="class",keep=TRUE,intercept=intc)
      
      alpha_mse[alphaID] = min(cv_glm$cvm)
      
      #cat('alphaID',alphaID,best_alpha_mse,'min',min(cv_glm$cvm),'\n')
      
      if(best_alpha_mse >= min(cv_glm$cvm)){
        #cat('changed\n')
        best_alpha = alphax
        best_alpha_mse = min(cv_glm$cvm)
        
      }
      
    }
    
  }
  
  cv_repeat = matrix(numeric(1),repeatRun,2)
  cv_permut = matrix(numeric(1),repeatRun,2)
  
  colnames(cv_repeat) = c('min','se1')
  colnames(cv_permut) = c('min','se1')
  
  for(testIdx in seq(1,repeatRun)){
    
    set.seed(testIdx)
    foldID = createFolds(func,k=cvFolds,list = F)
    
    cv_glm = cv.glmnet(feats,func,nfolds=cvFolds,foldid=foldID,family="multinomial",standardize=FALSE,weights=weights,alpha=best_alpha, type.multinomial = grpType,type.measure="class",keep=TRUE,intercept=intc)
    
    cv_perm = cv.glmnet(feats[sample(nrow(feats)),],func,nfolds=cvFolds,family="multinomial",standardize=FALSE,weights=weights,alpha=best_alpha, type.multinomial = grpType,type.measure="class",keep=TRUE,intercept=intc)
    
    cv_repeat[testIdx,'min'] = cv_glm$cvm[match(cv_glm$lambda.min,cv_glm$lambda)]
    cv_repeat[testIdx,'se1'] = cv_glm$cvm[match(cv_glm$lambda.1se,cv_glm$lambda)]
    
    cv_permut[testIdx,'min'] = cv_perm$cvm[match(cv_perm$lambda.min,cv_perm$lambda)]
    cv_permut[testIdx,'se1'] = cv_perm$cvm[match(cv_perm$lambda.1se,cv_perm$lambda)]
    
  }
  
  # Final run for visualilzation
  set.seed(8357)
  foldID = createFolds(func,cvFolds,list=F)
  
  cv_glm = cv.glmnet(feats,func,nfolds=cvFolds,foldid=foldID,family="multinomial",standardize=FALSE,weights=weights,alpha=best_alpha, type.multinomial = grpType,type.measure="class",keep=TRUE,intercept=intc)
  
  #save the prevalidated array for viz
  fit_preval_min = cv_glm$fit.preval[,,match(cv_glm$lambda.min,cv_glm$lambda)]
  fit_preval_1se = cv_glm$fit.preval[,,match(cv_glm$lambda.1se,cv_glm$lambda)]
  
  #get feature coefficients for the model
  cv_glm_coef_min = coef(cv_glm,s="lambda.min")
  cv_glm_coef_1se = coef(cv_glm,s="lambda.1se")
  
  mse_min = cv_glm$cvm[match(cv_glm$lambda.min,cv_glm$lambda)]
  mse_1se = cv_glm$cvm[match(cv_glm$lambda.1se,cv_glm$lambda)]
  
  numClasses = length(unique(func))
  
  betas_min = matrix(0,nrow=numFeat,ncol=numClasses)
  rownames(betas_min) = colnames(feats)
  betas_1se = matrix(0,nrow=numFeat,ncol=numClasses)
  rownames(betas_1se) = colnames(feats)
  
  for(classIdx in 1:numClasses){
    
    betas_min[,classIdx] = cv_glm_coef_min[[classIdx]][-1]
    betas_1se[,classIdx] = cv_glm_coef_1se[[classIdx]][-1]
    
  }
  
  min_model = list("lambda"=cv_glm$lambda.min,"mse"=mse_min,"preval"=fit_preval_min,"coeff"=betas_min)
  se1_model = list("lambda"=cv_glm$lambda.1se,"mse"=mse_1se,"preval"=fit_preval_1se,"coeff"=betas_1se)
  
  best_model = min_model
  
  #return(list("removed"=rmList,"final_func"=func,"alpha_mse"=alpha_mse,"best_alpha"=best_alpha,"final_fit"=cv_glm,"best_model"=best_model,"min_model"=min_model,"se1_model"=se1_model,"repeat_mse_min"=mean(cv_repeat[,'min']),"repeat_mse_se"=sd(cv_repeat[,'min']),"repeat_mse_se1"=mean(cv_repeat[,'se1']),"repeat_mse_sd"=sd(cv_repeat[,'se1']),"permut_mse_min"=mean(cv_permut[,'min']),"permut_mse_se"=sd(cv_permut[,'min']),"permut_mse_se1"=mean(cv_permut[,'se1']),"permut_mse_sd"=sd(cv_permut[,'se1'])))
  
  return(list("removed"=rmList,"final_func"=func,"alpha_mse"=alpha_mse,"best_alpha"=best_alpha,"final_fit"=cv_glm,"best_model"=best_model,"min_model"=min_model,"se1_model"=se1_model,"repeat_mse_min"=mean(cv_repeat[,'min']),"repeat_mse_min_sd"=sd(cv_repeat[,'min']),"repeat_mse_se1"=mean(cv_repeat[,'se1']),"repeat_mse_se1_sd"=sd(cv_repeat[,'se1']),"cv_repeat"=cv_repeat,"permut_mse_min"=mean(cv_permut[,'min']),"permut_mse_min_sd"=sd(cv_permut[,'min']),"permut_mse_se1"=mean(cv_permut[,'se1']),"permut_mse_se1_sd"=sd(cv_permut[,'se1']),"cv_permut"=cv_permut))
  
}