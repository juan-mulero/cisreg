#Function to perform sequence mapping between different reference genomes in a single step 
#through the call to the sequence of functions that must be executed to perform this task.
#Input: bed_table --> table in bed format (chr start end coord) with the original sequences to be remapped to another reference genome.
#       inputHg --> original reference genome. Format: integer 19 or 38
#       outputHg --> new reference genome. Format: integer 19 or 38
#       acceptDiffRatio --> acceptDiffRatio --> accepted size variation ratio between sequences, typically 0.05
hgtohg_onestep = function(bed_table, inputHg, outputHg, acceptDiffRatio){
  source("./scripts/hgtohg.R")
  table_hg2hg = hgtohg(bed_table, inputHg, outputHg)
  
  source("./scripts/size_comparison.R")
  size_hg2hg = size_comparison(bed_table, table_hg2hg, acceptDiffRatio)
  
  source("./scripts/filtAfterSizeComp.R")
  filtered = filtAfterSizeComp(table_hg2hg, size_hg2hg, bed_table)
  output = list(conversion = table_hg2hg, size_comparison = size_hg2hg, filtered = filtered)
  output
}