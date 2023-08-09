#Version of filtAfterSizeComp.R with parallel computing

filtAfterSizeComp_parallel = function(result_hgtohg, result_size_comp, orig_bedTable){
  indexes = c(result_size_comp$indexes_seq_notAcceptDiff, result_size_comp$indexes_NAs)
  index_NA = which(is.na(indexes))
  if (length(index_NA) > 0){indexes = indexes - index_NA}
  new_bedTable = result_hgtohg$table_oldHgtoNewHg
  merge_table = result_hgtohg$table_oldHgtoNewHg
  IDs = new_bedTable[c(4,8)]
  if (length(indexes) > 0){
    orig_bedTable = orig_bedTable[-indexes,]
    new_bedTable = new_bedTable[-indexes,]
    merge_table[indexes, 5:8] = NA
  }
  output = list(orig_bedTable = orig_bedTable[1:4], new_bedTable = new_bedTable[1:4], IDs = IDs, merge_table = merge_table)
}