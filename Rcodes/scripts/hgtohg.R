#Function to map sequences from one reference genome to another.
#Input: bed_table --> table in bed format (chr start end coord) with the original sequences to be remapped to another reference genome.
#       inputHg --> original reference genome. Format: integer 19 or 38
#       outputHg --> new reference genome. Format: integer 19 or 38
#Output: list with different objects:
#         table_oldHgtoNewHg: table with the sequences mapped to the new reference genome. 
#                            Also it includes the identifier (joined coordinates) of the original sequence and 
#                            the length of the segments that compose the sequence.
#         length_enh: vector with the number of segments that compose each sequence.
#         time: time of execution
hgtohg = function(bed_table, inputHg, outputHg){
  
  t = proc.time()
  
  #Packages
  library(liftOver)
  library(GenomicRanges)
  
  #Object GenomicRanges and LiftOver
  colnames(bed_table) = c("chr", "start", "end", "ID")
  over.chain = paste0("hg",inputHg,"ToHg",outputHg, ".over.chain")
  path = system.file(package = "liftOver", "extdata", over.chain)
  chain = import.chain(path)
  enhancers_liftover = GRanges(seqnames = bed_table$chr, 
                      ranges = IRanges(start = as.numeric(bed_table$start), end = as.numeric(bed_table$end)))
  enhancers_liftover = liftOver(enhancers_liftover, chain)
  
  #hg to hg
  n = 10000
  n_groups = ceiling(nrow(bed_table)/n)
  final_table_oldHgtoNewHg = data.frame()
  final_length_enh = c()
  init = 1
  for (group in 1:n_groups){
    if (group == n_groups){
      final = nrow(bed_table)
    } else {
      final = init + n - 1
    }
    
    #Initialization
    old_ID = new_ID = length_enh = chr = start = end = width_frag = c()
    
    #Table generator
    index = 1
    for (i in seq(init,final)){
      old_ID[index] = as.character(bed_table$ID[i])
      n_enhancers = length(enhancers_liftover[[i]])
      length_enh[index] = n_enhancers
      seq = enhancers_liftover[[i]]
      n_chr = length(seq@seqnames@values)
      
      if (n_enhancers == 1){
        new_enhancer = as.character(seq)
        new_ID[index] = new_enhancer
        chr[index] = as.character(seq@seqnames)
        start[index] = seq@ranges@start
        end[index] = seq@ranges@start + seq@ranges@width - 1
        width_frag[index] = seq@ranges@width
      } 
      else if (n_enhancers > 1 & n_chr == 1){
        start_enh = seq@ranges@start
        start_enh = min(start_enh)
        start[index] = start_enh
        
        end_enh = na.omit(seq@ranges@start + seq@ranges@width - 1)
        end_enh = max(end_enh) #NA values are omitted because they will be selected with de function max
        end[index] = end_enh
        
        chr_enh = unique(as.character(seq@seqnames@values))
        chr[index] = chr_enh
        
        width_frag[index] = sum(seq@ranges@width)
        
        new_enhancer = paste0(chr_enh, ":", start_enh, "-", end_enh)
        new_ID[index] = new_enhancer
      } 
      else {
        new_ID[index] = NA
        chr[index] = NA
        start[index] = NA
        end[index] = NA
        width_frag[index] = NA
      }
      index = index + 1
    }

    #Results
    table_oldHgtoNewHg = data.frame(chr, start, end, new_ID, old_ID, width_frag)
    final_table_oldHgtoNewHg = rbind(final_table_oldHgtoNewHg, table_oldHgtoNewHg)
    final_length_enh = c(final_length_enh, length_enh)
    init = final + 1
  }
  
  #Output
  t = proc.time() - t
  results = list(table_oldHgtoNewHg = final_table_oldHgtoNewHg, length_enh = final_length_enh, time = t)
  results
}

