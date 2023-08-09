#function to transform the crm coordinates of the tsv into the IDs used in the RDF. 
#The prior generation of IDs speeds up the execution of the cisreg go code.
#old_IDs --> vector with the hg38 coordinatates to convert (format: assembly_chr_start_end --> Ex: hg38_chr1_1234_5678)
#n --> splitting the identifier vector into groups of size n reduces execution time by reducing memory usage.

cisreg_IDs = function(old_IDs, n) {
  library(data.table)
  crm_new2old_IDs = fread("./crm/crm_new2old_IDs.tsv", header = F)
  colnames(crm_new2old_IDs) = c("new", "old")
  crm_old2new_IDs = fread("./crm/crm_old2new_IDs.tsv", header = F)
  colnames(crm_old2new_IDs) = c("old", "new")
  
  old_IDs = unique(old_IDs)
  old_IDs = gsub("hg|chr", "", old_IDs)
  old_IDs = gsub(":|-|_", ".", old_IDs)

  indexes = which(old_IDs %in% crm_old2new_IDs$old)
  if (length(indexes) > 0) {old_IDs = old_IDs[-indexes]}
  
  n_groups = ceiling(length(old_IDs)/n)
  init = 1
  for (group in 1:n_groups){
    if (group == n_groups){
      end = length(old_IDs)
    } else {
      end = init + n - 1
    }
      
    position = nrow(crm_new2old_IDs)
    new_IDs = c()
    for (i in init:end){
      position = position + 1
      n_zeros = 11 - nchar(position)
      new_ID = paste0("CRMHS", paste0(rep(0, n_zeros), collapse = ""), position)
      new_IDs = c(new_IDs, new_ID)
    }
      
    set = data.frame(old_IDs[init:end], new_IDs)
    set = format(set, scientific = F)
    colnames(set) = c("old", "new")
    crm_old2new_IDs = as.data.frame(crm_old2new_IDs)
    crm_old2new_IDs = rbind(crm_old2new_IDs, set)
      
    set = data.frame(new_IDs, old_IDs[init:end])
    colnames(set) = c("new", "old")
    crm_new2old_IDs = as.data.frame(crm_new2old_IDs)
    crm_new2old_IDs = rbind(crm_new2old_IDs, set)
      
    init = end + 1
  }

  write.table(crm_old2new_IDs, "./crm/crm_old2new_IDs.tsv", col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(crm_new2old_IDs, "./crm/crm_new2old_IDs.tsv", col.names = F, row.names = F, quote = F, sep = "\t")
  
  return("Done")
}
