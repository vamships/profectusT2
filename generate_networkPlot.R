# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : Profectus T2
# Description : Performs survival analysis using the biophysical measurements
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

dir_res = paste('results_networkPlot/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

num_folds = 8 # number of folds for cross-validation
featSelectMethod = 'upf' # Univariate polyserial filter
thresh_upf = 0.95 # polyserial coefficient threshold
stop_limit = 0.3 # log-likelihood cutoff for stopping backward search
num_repeat = 100 # number of repetitions of cross-validation
top_feat_thresh = 85 # Frequency cutoff to determine most-frequent features
seed_final = 8357 # seed to determine folds for final cross-validation model
doLookBack = T # To look for correlates of correlates after the final model
lookBackCorThresh = 0.75 # correlation cutoff to determine if features are correlated
set_plots = TRUE

log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

cat('\nCox-model parameters','\n',sep="",file=log_file,append=T)
cat('Number of folds : ',num_folds,'\n',file=log_file,append=T)
cat('Feature Selection Method : ',featSelectMethod,'\n',file=log_file,append=T)
cat('Threshold for upf : ',thresh_upf,'\n',file=log_file,append=T)
cat('Stop Limit for Backward Elimiation : ',stop_limit,'\n',file=log_file,append=T)
cat('Number of repeated evaluations : ',num_repeat,'\n',file=log_file,append=T)
cat('Threshold for frequency of features : ',top_feat_thresh,'\n',file=log_file,append=T)
cat('Seed for final model : ',seed_final,'\n',file=log_file,append=T)
cat('Activate lookBack feature : ',doLookBack,'\n',file=log_file,append=T)
cat('Correlation threshold for looking back : ',lookBackCorThresh,'\n',file=log_file,append=T)

dir_final = 'results_coxPred_biophysical/'

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

dir_fig_1 = paste(dir_res,'fig_1/',sep="")
dir.create(dir_fig_1)

dir_holistic = paste(dir_fig_1,'holistic/',sep="")
dir.create(dir_holistic)

feats_top = read.csv(paste(dir_final,'surv/final/feats.csv',sep=""),row.names = 1)
model_hol = coxSurvivalFinalModel(feats_top,subjects,seed_final,scolors,group_colors,km_colors,num_folds,set_plots,dir_holistic)

num_feat = ncol(feats_top)
num_cols = 7

nodes = as.data.frame(matrix(NA,nrow=num_feat,ncol=num_cols))
colnames(nodes) = c("name","reagent","antigen","reagent_color","antigen_color","cox_coeff","group")

coeffs_cox = round(log(model_hol$coeffs),2)

for(featIdx in 1:num_feat){
  
  featName = names(coeffs_cox)[featIdx]
  nodes[featIdx,"name"] = featName
  
  nodes[featIdx,"cox_coeff"] = coeffs_cox[featIdx]
  
  antigen_color = 'grey64'
  for(antigen in antigen_names){
    
    if(length(grep(antigen,featName))!=0){
      
      nodes[featIdx,"antigen"] = antigen
      nodes[featIdx,"antigen_color"] = antigen_colors[antigen]
      
    }
    
  }
  
  reagent_color = 'grey64'
  for(reagent in reagent_names){
    
    if(length(grep(reagent,featName))!=0){
      
      nodes[featIdx,"reagent"] = reagent
      nodes[featIdx,"reagent_color"] = reagent_colors[reagent]
      
    }
    
  }
  
  nodes[featIdx,"group"] = paste("g",featIdx,sep="")
  
}


# Substitution analysis
nodes_sub = NULL
feats_top = feats_top[,names(coeffs_cox)]

for(featIdx in 1:num_feat){
  
  featName = colnames(feats_top)[featIdx]
  feats_sel  = feats_top[,featIdx]
  xcor = cor(feats_sel,luminex,use='pairwise.complete.obs')
  
  otherFeatIdx = which(abs(xcor)>=0.75)
  
  if(length(otherFeatIdx)>1){
    
    feats_sub = feats_top[,-featIdx]
    colnames(feats_sub) = colnames(feats_top)[-featIdx]
    
    nodes_extra = as.data.frame(matrix(NA,nrow=length(otherFeatIdx),ncol=7))
    colnames(nodes_extra) = c("name","reagent","antigen","reagent_color","antigen_color","cox_coeff","group")
    
    for(subIdx in 1:length(otherFeatIdx)){
      
      subName = colnames(luminex)[otherFeatIdx[subIdx]]
      feats_add = cbind(feats_sub,luminex[,otherFeatIdx[subIdx]])
      colnames(feats_add) = c(colnames(feats_sub),subName)
      
      if(featName!=subName){
        
        model_hol = coxSurvivalFinalModel(feats_add,subjects,seed_final,scolors,group_colors,km_colors,num_folds,set_plots,dir_holistic)
        
        coeffs_cox = round(log(model_hol$coeffs),2)
        
        coeff_sub = coeffs_cox[subName]
        
        nodes_extra[subIdx,"name"] = subName
        nodes_extra[subIdx,"cox_coeff"] = coeff_sub
        
        antigen_color = 'grey64'
        for(antigen in antigen_names){
          
          if(length(grep(antigen,subName))!=0){
            
            nodes_extra[subIdx,"antigen"] = antigen
            nodes_extra[subIdx,"antigen_color"] = antigen_colors[antigen]
            
          }
          
        }
        
        reagent_color = 'grey64'
        for(reagent in reagent_names){
          
          if(length(grep(reagent,subName))!=0){
            
            nodes_extra[subIdx,"reagent"] = reagent
            nodes_extra[subIdx,"reagent_color"] = reagent_colors[reagent]
            
          }
          
        }
        
        nodes_extra[subIdx,"group"] = paste("g",featIdx,sep="")
        
      }
      
    }
    
    nodes_sub = rbind(nodes_sub,nodes_extra)
    
  }
  
}

rm_ID = which(is.na(nodes_sub[,"name"]))
if(length(rm_ID)!=0){
  
  nodes_sub = nodes_sub[-rm_ID,]
  
}

nodes = rbind(nodes,nodes_sub)

# make edges
num_feat = nrow(nodes)
num_rows = (num_feat * (num_feat-1))/2
edges = as.data.frame(matrix(NA,nrow=num_rows,ncol=6))
colnames(edges) = c("featName","nextName","source","target","pcc","pcc_pval")

edgeIdx = 0
for(featIdx in 1:(num_feat-1)){
  
  featName = nodes[featIdx,"name"]
  
  for(nextIdx in (featIdx+1):num_feat){
    
    nextName = nodes[nextIdx,"name"]
    
    edgeIdx = edgeIdx + 1
    edges[edgeIdx,"featName"] = featName
    edges[edgeIdx,"nextName"] = nextName
    edges[edgeIdx,"source"] = featIdx
    edges[edgeIdx,"target"] = nextIdx
    
    
    pcor_test = cor.test(luminex[,featName],luminex[,nextName])
    
    edges[edgeIdx,"pcc"] = pcor_test$estimate
    edges[edgeIdx,"pcc_pval"] = pcor_test$p.value
    
    
    
  }
  
}

file_out = paste(dir_res,"nodes.xml",sep="")
file.create(file_out)

cat("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n",file=file_out,append=T)
cat("<graph id=\"52\" label=\"Features\" directed=\"0\" cy:documentVersion=\"3.0\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:cy=\"http://www.cytoscape.org\" xmlns=\"http://www.cs.rpi.edu/XGMML\">\n",file=file_out,append=T)

for(nodeIdx in 1:nrow(nodes)){
  
  node_shape = "ELLIPSE"
  # if(nodes[nodeIdx,"cox_coeff"]>0){
  #   
  #   node_shape = "TRIANGLE"
  #   
  # }
  cat(paste("\n<node id=\"",nodeIdx,"\" label=\"",nodes[nodeIdx,"name"],"\">\n",sep=""),file=file_out,append=T)
  cat(paste("<att name=\"reagent\" value=\"",nodes[nodeIdx,"reagent"],"\" type=\"string\" cy:type=\"String\"/>\n",sep=""),file=file_out,append=T)
  cat(paste("<att name=\"antigen\" value=\"",nodes[nodeIdx,"antigen"],"\" type=\"string\" cy:type=\"String\"/>\n",sep=""),file=file_out,append=T)
  cat(paste("<att name=\"group\" value=\"",nodes[nodeIdx,"group"],"\" type=\"string\" cy:type=\"String\"/>\n",sep=""),file=file_out,append=T)
  cat(paste("<graphics h=\"65.0\" type=\"",node_shape,"\" width=\"20.0\" y=\"90.0\" outline=\"",col2hex(nodes[nodeIdx,"reagent_color"]),"\" w=\"65.0\" fill=\"",col2hex(nodes[nodeIdx,"antigen_color"]),"\">\n",sep=""),file=file_out,append=T)
  cat(paste("<att name=\"NODE_LABEL_POSITION\" value=\"","E,W,c,18.00,0.00","\" type=\"string\" cy:type=\"String\"/>\n",sep=""),file=file_out,append=T)
  cat(paste("<att name=\"NODE_LABEL_TRANSPARENCY\" value=\"","0","\" type=\"string\" cy:type=\"String\"/>\n",sep=""),file=file_out,append=T)
  cat(paste("</graphics>\n",sep=""),file=file_out,append=T)
  
  cat("</node>\n\n",file=file_out,append=T)
  
}

for(edgeIdx in 1:nrow(edges)){
  
  edge_color = "#00008B"
  if(edges[edgeIdx,"pcc"]<0){
    
    edge_color = "#B22222"
    
  }
  
  edge_type = "LONG_DASH"
  edge_transparency = "0"
  edge_width = "2"
  
  if(edges[edgeIdx,"pcc"]>0.75){
    
    edge_type = "SOLID"
    edge_width = "6"
    
    if(edges[edgeIdx,"pcc_pval"]<0.05){
      
      edge_transparency = "255"
      
    }
    
  }
  
  
  
  
  cat(paste("\n<edge id=\"",nodeIdx+edgeIdx,"\" label=\"",edges[edgeIdx,"featName"]," - ",edges[edgeIdx,"nextName"],"\" source=\"",edges[edgeIdx,"source"],"\" target=\"",edges[edgeIdx,"target"],"\">\n",sep=""),file=file_out,append=T)
  cat(paste("<graphics fill=\"",edge_color,"\" width=\"",edge_width,"\">\n",sep=""),file=file_out,append=T)
  cat(paste("<att name=\"EDGE_LINE_TYPE\" value=\"",edge_type,"\" type=\"string\" cy:type=\"String\"/>\n",sep=""),file=file_out,append=T)
  cat(paste("<att name=\"EDGE_TRANSPARENCY\" value=\"",edge_transparency,"\" type=\"string\" cy:type=\"String\"/>\n",sep=""),file=file_out,append=T)
  cat("</graphics>\n",file=file_out,append=T)
  
  cat("</edge>\n\n",file=file_out,append=T)
  
}

cat("</graph>\n",file=file_out,append=T)
