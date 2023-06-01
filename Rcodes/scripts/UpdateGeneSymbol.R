#Function to update gene symbols. For this task, first we select unique values to avoid duplicates, eliminate null values and 
#the genes with strange characters. Then we load Ensembl and run a query to get the updated gene symbols.
#Finally, we select those old gene symbols that ensembl has not been able to update and we make REST queries on the HGNC resource
#to obtain the updated gene. For this we search in aliases and previous symbols. 
#If a previous gene symbol currently corresponds to several genes, it is discarded, 
#because we will not discriminate which one corresponds to our particular case.
#Input: external_gene_name --> vector with the previous gene symbol
#Output: genes --> table with the previous gene symbol and the updated gene symbol
UpdateGeneSymbol = function(external_gene_name){
  library(biomaRt)
  source("./scripts/queryHGNC.R")
  
  #Pre-process: 
  external_gene_name = as.character(unique(external_gene_name))
  index = which(is.na(external_gene_name))
  if (length(index) != 0){external_gene_name = external_gene_name[-index]}
  external_gene_name = external_gene_name[grepl("^[[:alnum:][:blank:][:punct:]]+$", external_gene_name)]
  
  #Load Ensembl
  ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  #Query
  genes = getBM(attributes = c("external_gene_name", "hgnc_symbol"), filters = "external_gene_name",
                values = external_gene_name, mart = ensembl)
  genes[genes == ""] = NA
  genes = merge(genes, data.frame(external_gene_name), by = "external_gene_name", all = T)
  genes = genes[!duplicated(genes),]
  
  #Missing hgnc
  indexes_miss_hgnc = which(is.na(genes$hgnc_symbol))
  if (length(indexes_miss_hgnc) > 0){
    miss_hgnc = genes$external_gene_name[indexes_miss_hgnc]
    upd_miss_hgnc = queryHGNC(miss_hgnc)
    genes$hgnc_symbol[indexes_miss_hgnc] = upd_miss_hgnc
  }
  
  #Output
  genes
}