size_comparison_parallel = function(result_hgtohg, acceptDiffRatio){
  colnames(result_hgtohg$table_oldHgtoNewHg) = c("orig_chr", "orig_start", "orig_end", "orig_coord", 
                                                 "new_chr", "new_start", "new_end", "new_coord")
  old_size = as.numeric(as.character(result_hgtohg$table_oldHgtoNewHg$orig_end)) - as.numeric(
    as.character(result_hgtohg$table_oldHgtoNewHg$orig_start)) + 1
  new_size = as.numeric(as.character(result_hgtohg$table_oldHgtoNewHg$new_end)) - 
    as.numeric(as.character(result_hgtohg$table_oldHgtoNewHg$new_start)) + 1
  
  score = new_size/old_size
  indexes_NAs = which(is.na(new_size))
  NAs = result_hgtohg$table_oldHgtoNewHg$orig_coord[indexes_NAs]
  
  indexes_seq_diff_sizes = which(old_size != new_size)
  if (length(indexes_seq_diff_sizes) > 0) {
    size_diff = abs(new_size[indexes_seq_diff_sizes] - old_size[indexes_seq_diff_sizes])
    
    sizes_percent = c()
    for (i in 1:length(indexes_seq_diff_sizes)){
      index = indexes_seq_diff_sizes[i]
      accepted_size = acceptDiffRatio * old_size[index]
      if (size_diff[i] > accepted_size){sizes_percent[i] = F}
      else{sizes_percent[i] = T}
    }
    indexes = which(sizes_percent == F)
    indexes_seq_notAcceptDiff = indexes_seq_diff_sizes[indexes]
    
    seq_notAcceptDiff = data.frame(result_hgtohg$table_oldHgtoNewHg$orig_coord[indexes_seq_notAcceptDiff], 
                                   old_size[indexes_seq_notAcceptDiff], 
                                   result_hgtohg$table_oldHgtoNewHg$new_coord[indexes_seq_notAcceptDiff], 
                                   new_size[indexes_seq_notAcceptDiff], size_diff[indexes])
    colnames(seq_notAcceptDiff) = c("orig_coord", "old_size", "new_coord", "new_size", "size_diff")
    
    results = list(seq_notAcceptDiff = seq_notAcceptDiff, indexes_seq_notAcceptDiff = indexes_seq_notAcceptDiff, 
                   NAs = NAs, indexes_NAs = indexes_NAs, size_diff = size_diff, 
                   indexes_seq_diff_sizes = indexes_seq_diff_sizes, score = score)
  }
  else {
    results = list(seq_notAcceptDiff = data.frame(), indexes_seq_notAcceptDiff = c(), NAs = NAs, indexes_NAs = indexes_NAs, 
                   size_diff = c(), indexes_seq_diff_sizes = c(), score = score)
  }
  results
}