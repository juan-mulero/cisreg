#function to transform the crm coordinates of the tsv into the IDs used in the RDF. 
#The prior generation of IDs speeds up the execution of the cisreg go code.
#old_IDs --> vector with the hg38 coordinatates to convert (format: assembly_chr_start_end --> Ex: hg38_chr1_1234_5678)
#n --> splitting the identifier vector into groups of size n reduces execution time by reducing memory usage.

cisreg_IDs = function(old_IDs, n) {
  library(data.table)
  crm_old2new_IDs = fread("./crm/crm_old2new_IDs.tsv", header = F)
  colnames(crm_old2new_IDs) = c("old", "new")
  
  old_IDs = unique(old_IDs)
  old_IDs = gsub("hg|chr", "", old_IDs)
  old_IDs = gsub(":|-|_", ".", old_IDs)

  indexes = which(old_IDs %in% crm_old2new_IDs$old)
  if (length(indexes) > 0) {old_IDs = old_IDs[-indexes]}
  
  n_groups = ceiling(length(old_IDs)/n)
  init = 1
  initial_position = nrow(crm_old2new_IDs)
  rm(crm_old2new_IDs)
  new_crm_old2new_IDs = new_crm_new2old_IDs = c()
  t = proc.time()
  for (group in 1:n_groups){
    cat("Group", group, "of", n_groups, "with size of", n, "\n", proc.time() - t, "\n")
    if (group == n_groups){
      end = length(old_IDs)
    } else {
      end = init + n - 1
    }
    
    if (is.null(nrow(new_crm_old2new_IDs))){
      position = initial_position
    } else {position = initial_position + nrow(new_crm_old2new_IDs)}

    new_IDs = c()
    for (i in init:end){
      position = position + 1
      n_zeros = 11 - nchar(position)
      new_ID = paste0("CRMHS", paste0(rep(0, n_zeros), collapse = ""), position)
      new_IDs = c(new_IDs, new_ID)
    }
      
    set = data.table(old_IDs[init:end], new_IDs)
    set = format(set, scientific = F)
    colnames(set) = c("old", "new")
    new_crm_old2new_IDs = as.data.table(new_crm_old2new_IDs)
    new_crm_old2new_IDs = rbind(new_crm_old2new_IDs, set)
      
    set = data.table(new_IDs, old_IDs[init:end])
    colnames(set) = c("new", "old")
    new_crm_new2old_IDs = as.data.table(new_crm_new2old_IDs)
    new_crm_new2old_IDs = rbind(new_crm_new2old_IDs, set)
      
    init = end + 1
  }
  
  crm_old2new_IDs = fread("./crm/crm_old2new_IDs.tsv", header = F)
  colnames(crm_old2new_IDs) = colnames(new_crm_old2new_IDs)
  new_crm_old2new_IDs$old = as.character(new_crm_old2new_IDs$old)
  new_crm_old2new_IDs$new = as.character(new_crm_old2new_IDs$new)
  crm_old2new_IDs = rbind(crm_old2new_IDs, new_crm_old2new_IDs)
  write.table(crm_old2new_IDs, "./crm/crm_old2new_IDs.tsv", col.names = F, row.names = F, quote = F, sep = "\t")
  rm(crm_old2new_IDs, new_crm_old2new_IDs)
  
  crm_new2old_IDs = fread("./crm/crm_new2old_IDs.tsv", header = F)
  colnames(crm_new2old_IDs) = colnames(new_crm_new2old_IDs)
  new_crm_new2old_IDs$new = as.character(new_crm_new2old_IDs$new)
  new_crm_new2old_IDs$old = as.character(new_crm_new2old_IDs$old)
  crm_new2old_IDs = rbind(crm_new2old_IDs, new_crm_new2old_IDs)
  write.table(crm_new2old_IDs, "./crm/crm_new2old_IDs.tsv", col.names = F, row.names = F, quote = F, sep = "\t")
  rm(crm_new2old_IDs, new_crm_new2old_IDs)
  
  return("Done")
}
