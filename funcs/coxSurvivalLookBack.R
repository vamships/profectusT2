# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Decription  : Find correlates of the final features identified in survival analysis
# Cite        : TBD
# ******************************************************************************

coxSurvivalLookBack = function(feats,feats_all,top_feat_idx,lookBackCorThresh,baselineCindex,seedMe,subj_original,scolors_original,lcolors_all,group_colors,km_colors,num_folds,log_file,dir_lookBack){
  
  feats_top  = feats[,top_feat_idx]
  xcor = cor(feats_top,feats_all,use='pairwise.complete.obs')
  for(featIdx in 1:length(top_feat_idx)){
    
    # shortlist features with which the feature is correlated above the threshold
    otherFeats = which(abs(xcor[featIdx,])>=lookBackCorThresh)
    
    if(length(otherFeats)>1){
      
      subConcIdx = numeric(length(otherFeats))
      subFeatNames = character(length(otherFeats))
      subColors = character(length(otherFeats))
      
      dir_feat = paste(dir_lookBack,featIdx,'.',colnames(feats_top)[featIdx],"/",sep="")
      dir.create(dir_feat)
      
      feats_sub = feats_top[,-featIdx]
      
      # Perform cross-validation by substituting the original feature with each shortlisted feature
      otherFeatCount = 0
      for(otherFeatIdx in otherFeats){
        
        otherFeatCount = otherFeatCount + 1
        subFeatNames[otherFeatCount] = colnames(feats_all)[otherFeatIdx]
        subColors[otherFeatCount] = lcolors_all[otherFeatIdx,1]
        
        dir_otherFeat = paste(dir_feat,otherFeatIdx,'.',colnames(feats_all)[otherFeatIdx],"/",sep="")
        dir.create(dir_otherFeat)
        
        cat('replacing ',colnames(feats_top)[featIdx],' with ',colnames(feats_all)[otherFeatIdx],'\n')
        
        feats_add = data.frame(cbind(feats_sub,feats_all[,otherFeatIdx]))
        colnames(feats_add) = c(colnames(feats_sub),colnames(feats_all)[otherFeatIdx])
        
        write.csv(feats_add,file=paste(dir_otherFeat,'feats.csv',sep=""),row.names = T,col.names = T)
        
        model_add = coxSurvivalFinalModel(feats_add,subj_original,seedMe,scolors_original,group_colors,km_colors,num_folds,plots=TRUE,dir_otherFeat)
        
        subConcIdx[otherFeatCount] = model_add$cindex_test
        
        cat('\nreplacing ',colnames(feats_top)[featIdx],' with ',colnames(feats_all)[otherFeatIdx],'\n',file=log_file,append=T)
        cat('Replaced : ',model_add$cindex_train,'(Train)',model_add$cindex_test,'(Test)','\n',file=log_file,append=T)
        
      }
      
      pdf(paste(dir_feat,'subst_',featIdx,'.pdf',sep=""))
      par(mar=c(12,4,2,0.5))
      bp = barplot(subConcIdx,main=paste("Concordance Index for substituting\n",colnames(feats_top)[featIdx],sep=""),ylim=c(0,1),las=2,names.arg=subFeatNames,cex.names=0.6,col=subColors)
      abline(h=baselineCindex,lty=2,lwd=2)
      text(x=bp,y=subConcIdx-0.2,labels=subConcIdx,cex=0.8,xpd=TRUE)
      dev.off()
      
    }
    
  }
  
}