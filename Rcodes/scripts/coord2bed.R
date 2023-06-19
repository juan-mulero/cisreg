#Function to obtain a table in bed format from genomic coordinates. 
#The genomic coordinates follow the format chrX:start-end
#Input: coord --> vector with strings with the commented format
#Output: bedTable --> data.frame with bed format (chr  start end  coord)
coord2bed = function(coord){
  chr = start = end = c()
  for (i in 1:length(coord)){
    if (!is.na(coord[i])){
      split = unlist(strsplit(coord[i], "[:-]"))
      chr[i] = split[1]
      start[i] = as.numeric(split[2])
      end[i] = as.numeric(split[3])
    } else {
      chr[i] = NA
      start[i] = NA
      end[i] = NA
    }
  }
  bed = data.frame(chr, start, end, coord)
  bed
}