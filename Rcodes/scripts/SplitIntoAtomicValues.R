#Function to obtain a table with atomic values when the original table has a column with multiple values separated by a 
#specific character. The function splits this variable into its atomic values and generates individual lines for each value.
#Input:   table --> table to split
#         column2split --> number of the column to split
#         marksplit --> mark that split the different values.
#Output:  output -- > table with the selected column split in atomic values.
#Note: This function is only useful if the values of the other variables are linked equally to each of these multiple values. 
#Do not use if another variable also has multiple values and there is a direct relationship between 
#the position of these multiple values.
SplitIntoAtomicValues = function(table, column2split, marksplit){
  table = data.frame(table)
  colname = colnames(table)[column2split]
  indexes = grep(marksplit, table[,column2split])
  if (length(indexes) > 0){
    set1 = table[-indexes,]
    set2 = table[indexes,]
    row_ID = seq(1:nrow(set2))
    
    column = set2[[column2split]]
    new_row_ID = new_column = c()
    for (i in 1:length(column)){
      split = unlist(strsplit(column[i], marksplit))
      new_column = c(new_column, split)
      IDs = rep(row_ID[i], length(split))
      new_row_ID = c(new_row_ID, IDs)
    }
    
    new_table = cbind(row_ID = new_row_ID, new_column)
    new_table = new_table[!duplicated(new_table),]
    set2 = cbind(row_ID, set2[,-column2split])
    set2 = merge(set2, new_table, by = "row_ID", all.x = T)
    set2 = set2[!duplicated(set2),]
    set2 = set2[,-1]
    if (ncol(set2) == column2split) {
      set2 = set2[,c(1:(column2split-1), ncol(set2))]
    } else if (column2split == 1) {
      set2 = set2[,c(ncol(set2), 1:(ncol(set2)-1))]
    } else {set2 = set2[,c(1:(column2split-1), ncol(set2), column2split:(ncol(set2)-1))]}
    colnames(set2)[column2split] = colname
    table = rbind(set1, set2)
  }
  return(table)
}
