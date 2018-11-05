# color stuff (modified from Karen's Mario code)

# createColumnColors = function(featNames,reagent_colors,antigen_colors,strain_colors){
#   
#   # create column colors
#   lcolors = matrix(nrow=length(featNames), ncol=3, dimnames=list(featNames, c('reagent','antigen','strain')))
#   
#   lcolors[grep("IgG.",featNames),'reagent'] = reagent_colors[1]
#   lcolors[grep("R3A.",featNames),'reagent'] = reagent_colors[2]
#   lcolors[grep("R2A.",featNames),'reagent'] = reagent_colors[3]
#   lcolors[grep("MBL.",featNames),'reagent'] = reagent_colors[4]
#   lcolors[grep("FcgRIIa.",featNames),'reagent']= reagent_colors[5]
#   lcolors[grep("FcgRIIIa.",featNames),'reagent'] = reagent_colors[6]
#   lcolors[grep("FcgRIIIb.",featNames),'reagent'] = reagent_colors[7]
#   lcolors[grep("C1q",featNames),'reagent'] = reagent_colors[8]
#   
#   lcolors[grep("gp120",featNames),'antigen'] = antigen_colors[1]
#   lcolors[grep("gp130",featNames),'antigen'] = antigen_colors[2]
#   lcolors[grep("gp140",featNames),'antigen'] = antigen_colors[3]
#   lcolors[grep("gp160",featNames),'antigen'] = antigen_colors[4]
#   lcolors[grep("PR55",featNames),'antigen'] = antigen_colors[5]
#   lcolors[grep("cv2a",featNames),'antigen'] = antigen_colors[6]
#   lcolors[grep("cv2c",featNames),'antigen'] = antigen_colors[7]
#   
#   lcolors[grep("Ertl",featNames),'strain'] = strain_colors[1]
#   lcolors[grep("SIVcpz",featNames),'strain'] = strain_colors[2]
#   lcolors[grep("Novartis",featNames),'strain'] = strain_colors[3]
#   lcolors[grep("SIVmac1A11",featNames),'strain'] = strain_colors[4]
#   lcolors[grep("SIVmac239",featNames),'strain'] = strain_colors[5]
#   lcolors[grep("SIVsmE543",featNames),'strain'] = strain_colors[6]
#   lcolors[grep("SIVsmH4",featNames),'strain'] = strain_colors[7]
#   
#   return(lcolors)
#   
# }

createColumnColors = function(featNames,reagent_names,reagent_colors,antigen_names,antigen_colors){
  
  # create column colors
  lcolors = matrix(nrow=length(featNames), ncol=2, dimnames=list(featNames, c('reagent','antigen')))
  
  for(reagent_name in reagent_names){
    
    lcolors[grep(reagent_name,featNames),'reagent'] = reagent_colors[reagent_name]
    
  }
  
  for(antigen_name in antigen_names){
    
    lcolors[grep(antigen_name,featNames),'antigen'] = antigen_colors[antigen_name]
    
  }
  
  return(lcolors)
  
}