# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectsu T2
# Decription  : Helper function for backward feature selection for survival analysis
# Cite        : TBD
# ******************************************************************************

doCoxBackSearch = function(feats,survObject,selected_feat_idx,stop_limit,plots,dir_res){
  
  #cat('In doCoxBackSearch\n')
  
  feats_temp = feats
  feat_idx_list = selected_feat_idx
  
  path = vector(mode='numeric',length=ncol(feats))
  
  pathIdx = 0
  stop_search = FALSE
  
  current_model = coxph(survObject ~ .,data=feats_temp,control=coxph.control(iter.max=100))
  
  LL_prev = current_model$loglik[2]
  
  while(!stop_search & (ncol(feats_temp)>2)){
    
    current_model = coxph(survObject~.,data=feats_temp,control=coxph.control(iter.max=100))
    
    sel_sub = takeOffOneFeat(feats_temp,survObject,LL_prev,feat_idx_list,plots,dir_res)  #take off one feature
    
    pathIdx = pathIdx + 1
    path[pathIdx] = sel_sub$res['LL_opt']
    
    if(path[pathIdx]>=LL_prev | abs(path[pathIdx]-LL_prev) < round(log2(1/(1-stop_limit)),digits=2)){
      
      feats_temp = sel_sub$feats
      feat_idx_list = sel_sub$feat_idx_list
      
    }else{
      
      stop_search=TRUE
      
    }
    
  }
  
  path = path[path!=0]
  
  if(plots){
  
  pdf(paste(dir_res,'backSearch.pdf',sep=""))
  plot(-1:length(path),c(current_model$loglik[1],LL_prev,path),type='b',main='Backward Search Feature Selection',xlab='step #',ylab='Likelihood')
  dev.off()
  
  pdf(paste(dir_res,'backSearchtry.pdf',sep=""))
  xTickLabels = c('null','all',1:length(path))
  plot(-1:length(path),c(current_model$loglik[1],LL_prev,path),type='b',main='Backward Search Feature Selection',xlab='step #',ylab='Likelihood',xaxt='n')
  axis(1,labels=FALSE)
  text(-1:length(path),par('usr')[3]-0.4,srt=90,adj=1,labels=xTickLabels,xpd=TRUE)
  dev.off()
  
  }
  
  return(feat_idx_list)
  
}