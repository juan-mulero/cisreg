#Function to filter out those sequences that could not be remapped with another reference genome.
#This function must be executed after the previous mapping (hgtohg) and after the length comparison (size_comparison).
#Input: result_hgtohg --> result of the function hgtohg.
#       result_size_comp --> result of the function size_comparison.
#       orig_bedTable --> table in bed format with the original sequences to be remapped to another reference genome.
#Output: list with different objects:
#         orig_bedTable --> table in bed format with the original sequences to be remapped to another reference genome.
#         new_bedTable --> table in bed format with the new sequences remapped to another reference genome.
#         IDs --> table linking the coordinates in both reference genomes.
filtAfterSizeComp = function(result_hgtohg, result_size_comp, orig_bedTable){
  indexes = c(result_size_comp$indexes_seq_notAcceptDiff, result_size_comp$indexes_NAs)
  new_bedTable = result_hgtohg$table_oldHgtoNewHg
  IDs = new_bedTable[4:5]
  if (length(indexes) > 0){
    orig_bedTable = orig_bedTable[-indexes,]
    new_bedTable = new_bedTable[-indexes,]
  }
  colnames(orig_bedTable)[4] = "old_ID"
  output = list(orig_bedTable = orig_bedTable, new_bedTable = new_bedTable[1:4], IDs = IDs)
}
