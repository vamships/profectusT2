# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectus T2
# Decription  : Implements pre-filtering based on polyserial correlation with time-to-infection
# Cite        : TBD
# ******************************************************************************

univariatePolyserialFilter = function(corr_thresh,feats,subj_original,scolors,lcolors,plots,dir_slope){
  
  # Compute polyserial coefficient of each feature
  xcor = matrix(NA,nrow=ncol(feats),ncol=1)
  colnames(xcor) = c('Overall')
  
  for(featIdx in 1:ncol(feats)){
    
    xcor[featIdx,1] = polyserial(feats[,featIdx],subj_original[,'Challenges'])
    
  }
  
  if(plots){
    
    pdf(paste(dir_slope,'poly_dist_1.pdf',sep=""))
    yrange = range(xcor,na.rm=TRUE)
    plot(c(1,ncol(feats)),c(yrange[1],yrange[2]),type='n',main='Polyserial Coefficients',xlab='Features',ylab='')
    lines(1:ncol(feats),xcor[,1],type='l',lty=1,col='black',lwd=3)
    dev.off()
    
    pdf(paste(dir_slope,'poly_dist_2.pdf',sep=""))
    boxplot(abs(xcor),main='Polyserial Coefficients')
    dev.off()
    
    pdf(paste(dir_slope,'poly_dist_overall.pdf',sep=''))
    plot(1:nrow(xcor),abs(xcor[,1]),type='b',main='Polyserial coeffs - Overall',xlab='Feature index',ylab='Coeff. Value',col='black')
    abline(h=quantile(abs(xcor[,1]),corr_thresh),lty=2,col='black')
    legend('topleft',legend=paste('thresh : ',corr_thresh,sep=''),lty=2)
    dev.off()
    
  }
  
  
  sortFeats = function(pcor,selected_feat_idx){
    
    # View polyserial filter features
    pcor_select = pcor[selected_feat_idx]
    pcor_order = order(pcor_select,decreasing = TRUE,na.last=TRUE)
    selected_feat_idx = selected_feat_idx[pcor_order]
    
    return(selected_feat_idx)
    
  }
  
  
  #Select based on magnitude of polyserial coefficient
  selected_feat_idx_ovrl = which(abs(xcor[,1])>quantile(abs(xcor[,1]),corr_thresh))
  selected_feat_idx_ovrl = sortFeats(abs(xcor[,1]),selected_feat_idx_ovrl)
  
  cat('preliminary subset for overall- ',selected_feat_idx_ovrl,'\n')
  
  feats_sel = feats[,selected_feat_idx_ovrl]
  
  if(plots){
    
    pdf(paste(dir_slope,'selected_overall.pdf',sep=""))
    par(mar=c(12,4,2,0.5))
    barplot(xcor[selected_feat_idx_ovrl,1],main='Polyserials of selected features - Overall',ylab='Corr. Coeff',las=2,names.arg=colnames(feats[,selected_feat_idx_ovrl]),cex.names=0.7)
    dev.off()
    
    corr_feat_sel = cor(feats_sel,use='pairwise.complete.obs')
    pdf(paste(dir_slope,'int_full_corr_feat_all.pdf',sep=""))
    corrplot(corr_feat_sel,tl.cex = 0.7,addCoef.col = "black")
    dev.off()
    
  }
  
  #Reduce pair-wise correlations between selected features
  mark_sel = matrix(0,nrow=1,ncol=ncol(feats_sel))
  rownames(mark_sel) = 'selection'
  colnames(mark_sel) = colnames(feats_sel)
  
  for(featIdx in 1:ncol(feats_sel)){
    
    if(mark_sel[1,featIdx]==0){
      
      # Mark as selected
      mark_sel[1,featIdx] = 1
      
      # Feats 
      idx_to_consider = which(mark_sel==0)
      
      if(length(idx_to_consider)>0){
        
        cor_consider = cor(feats_sel[,featIdx],feats_sel[,idx_to_consider],use='pairwise.complete.obs')
        idx_throw = which(cor_consider>0.8)
        
        if(length(idx_throw)>0){
          
          mark_sel[1,idx_to_consider[idx_throw]] = 2
          
        }
        
      }
      
    }
    
  }
  
  sub_feat_idx_ovrl = which(mark_sel==1)
  
  # if selected less than one, select the next best also
  if(length(sub_feat_idx_ovrl)<2){
    
    sub_feat_idx_ovrl = c(sub_feat_idx_ovrl,2)
    
  }
  
  cat("Initial Selected:",length(selected_feat_idx_ovrl),"; After reducing pair-wise correlations:",length(sub_feat_idx_ovrl),"\n")
  
  if(plots){
    
    pdf(paste(dir_slope,'selected_overall.pdf',sep=""))
    par(mar=c(12,4,2,0.5))
    barplot(xcor[selected_feat_idx_ovrl,1],main='Polyserials of selected features - Overall',ylab='Corr. Coeff',las=2,names.arg=colnames(feats[,selected_feat_idx_ovrl]),cex.names=0.7,addCoef.col = "black")
    dev.off()
    
    if(length(sub_feat_idx_ovrl)>1){
      corr_sub_sel = cor(feats_sel[,sub_feat_idx_ovrl],use='pairwise.complete.obs')
      pdf(paste(dir_slope,'int_sel_corr_feat_all.pdf',sep=""))
      corrplot(corr_sub_sel,tl.cex=0.7,addCoef.col = "black")
      dev.off()
    }
    
    makeSimpleHeatMaps = function(feats_select,subj_select,scolors_select,lcolors_select,fname){
      
      chall.sort = sort(subj_select[,'Challenges'], index.return=TRUE)
      
      ldata = scale(feats_select)
      ldata[ldata>3] = 3; ldata[ldata < -3] = -3
      lr = ceiling(max(abs(min(ldata,na.rm=TRUE)), max(ldata,na.rm=TRUE)))
      lbreaks = seq(-lr,lr,0.1)
      
      pdf(paste(fname,'-raw.pdf'))
      par(mar=c(12,4,2,0.5))
      heatmap.4(ldata, col=bluered, scale='none', trace='none', cexRow=0.5, cexCol=0.8, margin=c(8,5), breaks=lbreaks,  key=FALSE,symkey=T, dendrogram='none',hclust=hclust.ward,RowSideColors=scolors_select,NumRowSideColors=4,ColSideColors=as.matrix(lcolors_select),NumColSideColors=2,Rowv=FALSE, Colv=FALSE, na.color='black')
      dev.off()
      
      pdf(paste(fname,'-by-challenge.pdf'))
      par(mar=c(12,4,2,0.5))
      heatmap.4(ldata[chall.sort$ix,], col=bluered, scale='none', trace='none', cexRow=0.5, cexCol=0.8, margin=c(8,5), breaks=lbreaks, key=T,symkey=FALSE, dendrogram='none',hclust=hclust.ward,RowSideColors=scolors_select[chall.sort$ix,],NumRowSideColors=4,ColSideColors=as.matrix(lcolors_select),NumColSideColors=2,Rowv=FALSE, Colv=FALSE, na.color='black')
      dev.off()
      
    }
    
    #
    lcolors_ovrl = as.matrix(colorRampPalette(c('goldenrod3','lightgoldenrod1'))(length(selected_feat_idx_ovrl)))
    lcolors_ovrl = cbind(lcolors_ovrl,lcolors_ovrl)
    lcolors_ovrl[sub_feat_idx_ovrl,1] = "green"
    makeSimpleHeatMaps(feats[,selected_feat_idx_ovrl],subj_original,scolors,lcolors_ovrl,paste(dir_slope,'luminex_overall',sep=""))
    
  }
  
  return(selected_feat_idx_ovrl[sub_feat_idx_ovrl])
  
}