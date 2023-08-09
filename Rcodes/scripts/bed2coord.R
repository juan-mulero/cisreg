#Function to obtain genomic coordinates from a table in bed format. 
#The genomic coordinates produced follow the format chrX:start-end
#Input: bedTable --> data.frame with bed format (chr  start end)
#Output: coord --> vector with strings with the commented format
bed2coord = function(bedTable){
  colnames(bedTable)[1:3] = c("chr", "start", "end")
  bedTable = dplyr::mutate(bedTable, coord = paste0(chr, ":", start, "-", end))
  coord = bedTable$coord
}