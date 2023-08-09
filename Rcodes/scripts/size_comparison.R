#Function to compare the length of the sequences after their mapping between reference genomes.
#Input: origBedTable --> table in bed format used to map the sequences to the new reference genome.
#       result_hgtog --> result of the function hgtohg.
#       acceptDiffRatio --> acceptDiffRatio --> accepted size variation ratio between sequences, typically 0.05
#Output: list with different objects:
#       seq_notAcceptDiff --> table with the sequences that do not satisfy the condition of the accepted variation ratio.
#       indexes_seq_notAcceptDiff --> indexes of these sequences in the original table of sequences.
#       NAs --> vector with the sequences that could not be mapped for reasons different to the size
#       indexes_NAs --> indexes of these sequences in the original table of sequences.
#       size_diff --> length difference in the sequences that satisfy the accepted variation ratio.
#       indexes_seq_diff_sizes --> indexes of these sequences in the original table of sequences.
#       score --> score between both sizes
size_comparison = function(origBedTable, result_hgtohg, acceptDiffRatio){
  colnames(origBedTable) = c("chr", "start", "end", "ID")
  old_size = as.numeric(as.character(origBedTable$end)) - as.numeric(as.character(origBedTable$start)) + 1
  new_size = as.numeric(as.character(result_hgtohg$table_oldHgtoNewHg$end)) - 
             as.numeric(as.character(result_hgtohg$table_oldHgtoNewHg$start)) + 1
  
  score = new_size/old_size
  indexes_NAs = which(is.na(new_size))
  NAs = result_hgtohg$table_oldHgtoNewHg$old_ID[indexes_NAs]
  
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
    
    seq_notAcceptDiff = data.frame(result_hgtohg$table_oldHgtoNewHg$old_ID[indexes_seq_notAcceptDiff], 
                                   old_size[indexes_seq_notAcceptDiff], 
                                   result_hgtohg$table_oldHgtoNewHg$new_ID[indexes_seq_notAcceptDiff], 
                                   new_size[indexes_seq_notAcceptDiff], size_diff[indexes])
    colnames(seq_notAcceptDiff) = c("old_ID", "old_size", "new_ID", "new_size", "size_diff")
    
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