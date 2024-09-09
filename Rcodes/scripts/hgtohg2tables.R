#Function to obtain tables with a format of interest to build the final file.
#Input: oldhg --> original assembly. Format: "hg19" or "hg38"
#       newhg --> current assemby. Formta: "hg19" or "hg38"
#       result_hgtohg_onestep --> result of the function hgtohg_onestep
#       table_oldhg --> table with the original sequences.
#Output: table_oldhg --> table with the original sequences.
#        table_newhg --> table with the new sequences.
#        enhancers --> table with the merge between both previous tables.
hgtohg2tables = function(oldhg, newhg, table_oldhg, result_hgtohg_onestep){
  table_newhg = result_hgtohg_onestep$filtered$new_bedTable
  colnames(table_newhg)[4] = "coord"
  table_newhg = cbind(table_newhg, current_assembly = newhg)
 
  enh_coord = result_hgtohg_onestep$filtered$IDs
  colnames(enh_coord) = c("coord", oldhg)
  enh_coord = cbind(enh_coord, minimum_ratio = 0.95, score = result_hgtohg_onestep$size_comparison$score)
  
  table_newhg = merge(table_newhg, enh_coord, by = "coord", all.x = T)
  table_newhg = table_newhg[!duplicated(table_newhg),]
  colnames(table_newhg)[c(1:4,6)] = c(newhg, "current_chr", "current_start", "current_end", "coord")
  table_newhg = merge(table_newhg, table_oldhg[,c(1,5)], by = "coord", all.x = T)
  table_newhg = table_newhg[,c(9,3:5,2,6:8)]
  colnames(table_newhg)[5] = "coord"
  
  enhancers = merge(table_oldhg[,-5], table_newhg[,-5], by = "enh_ID", all = T)
  enhancers = enhancers[!duplicated(enhancers),]
  
  output = list(table_oldhg = table_oldhg, table_newhg = table_newhg, enhancers = enhancers)
  return(output)
}
