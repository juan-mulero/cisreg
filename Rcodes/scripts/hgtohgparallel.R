#Version of hgtohg.R with parallel computing
#Function to map sequences from one reference genome to another.
#Input: bed_table --> table in bed format (chr start end ID) with the original sequences to be remapped to another reference genome.
#       inputHg --> original reference genome. Format: integer 19 or 38
#       outputHg --> new reference genome. Format: integer 19 or 38
#Output: list with different objects:
#         table_oldHgtoNewHg: table with the sequences mapped to the new reference genome. 
#                            Also it includes the identifier (joined coordinates) of the original sequence and 
#                            the length of the segments that compose the sequence.
#         time: time of execution
hgtohgparallel = function(bed_table, inputHg, outputHg, ncl){
  
  t = proc.time()

  #Packages
  library(liftOver)
  library(GenomicRanges)
  library(foreach)
  library(doParallel)
  
  #Object GenomicRanges and LiftOver
  colnames(bed_table) = c("chr", "start", "end", "ID")
  over.chain = paste0("hg",inputHg,"ToHg",outputHg, ".over.chain")
  path = system.file(package = "liftOver", "extdata", over.chain)
  chain = import.chain(path)
  enhancers_liftover = GRanges(seqnames = bed_table$chr, 
                               ranges = IRanges(start = as.numeric(bed_table$start), end = as.numeric(bed_table$end)))
  enhancers_liftover = liftOver(enhancers_liftover, chain)
  
  #Clusters
  cl = makeCluster(ncl)
  registerDoParallel(cl)
  
  #Table generator
  coord = NULL
  try({
    coord = foreach(i = 1:ncl, .combine = rbind) %dopar%{
      n_enhancers = length(enhancers_liftover[[i]])
      seq = enhancers_liftover[[i]]
      n_chr = length(seq@seqnames@values)
      
      if (n_enhancers == 1){new_ID = as.character(seq)
      } else if (n_enhancers > 1 & n_chr == 1){
        start_enh = seq@ranges@start
        start_enh = min(start_enh)
        
        end_enh = na.omit(seq@ranges@start + seq@ranges@width - 1)
        end_enh = max(end_enh) #NA values are omitted because they will be selected with de function max
        
        chr_enh = unique(as.character(seq@seqnames@values))
        
        new_ID = paste0(chr_enh, ":", start_enh, "-", end_enh)
      } else {new_ID = NA}
      
      new_ID
    }
  }, silent = T)
  
  #Table generator
  coord = foreach(i = 1:length(enhancers_liftover), .combine = rbind) %dopar%{
    n_enhancers = length(enhancers_liftover[[i]])
    seq = enhancers_liftover[[i]]
    n_chr = length(seq@seqnames@values)
    
    if (n_enhancers == 1){new_ID = as.character(seq)
    } else if (n_enhancers > 1 & n_chr == 1){
      start_enh = seq@ranges@start
      start_enh = min(start_enh)
      
      end_enh = na.omit(seq@ranges@start + seq@ranges@width - 1)
      end_enh = max(end_enh) #NA values are omitted because they will be selected with de function max
      
      chr_enh = unique(as.character(seq@seqnames@values))
      
      new_ID = paste0(chr_enh, ":", start_enh, "-", end_enh)
    } else {new_ID = NA}
    
    new_ID
  }
  
  #Close clusters
  stopCluster(cl)
  
  #Final table
  source("./scripts/coord2bed.R")
  new_bed = coord2bed(coord)
  colnames(new_bed) = c("new_chr", "new_start", "new_end", "new_coord")
  colnames(bed_table) = c("orig_chr", "orig_start", "orig_end", "orig_coord")
  table_oldHgtoNewHg = cbind(bed_table, new_bed)

  t = proc.time() - t
  
  output = list(table_oldHgtoNewHg = table_oldHgtoNewHg, time = t)
  output
}





