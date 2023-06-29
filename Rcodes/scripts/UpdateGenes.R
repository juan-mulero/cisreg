#Same as UpdateGeneSymbol.R, but performing queries on the .txt file instead of REST queries (faster)

UpdateGenes = function(external_gene_name){
  library(biomaRt)
  
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
    library(data.table)
    miss_hgnc = genes$external_gene_name[indexes_miss_hgnc]
    old_miss_hgnc = upd_miss_hgnc = c()
    
    temp = tempfile()
    download.file("http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt", temp)
    set_hgnc = fread(temp, header = T, quote = "")
    set_hgnc = data.frame(set_hgnc)
    rm(temp)
    
    for (j in 1:length(miss_hgnc)){
      patterns = c(paste0("^", miss_hgnc[j], "$"), paste0("^", miss_hgnc[j], "\\|"), 
                   paste0("\\|", miss_hgnc[j], "$"), paste0("\\|", miss_hgnc[j], "\\|"))
      patterns = paste(patterns, collapse = "|")
      symbol = grep(patterns, set_hgnc$symbol)
      prev_symbol = grep(patterns, set_hgnc$prev_symbol)
      alias_symbol = grep(patterns, set_hgnc$alias_symbol)
      
      if (length(symbol) > 0){
        upd_miss_hgnc = c(upd_miss_hgnc, set_hgnc$symbol[symbol])
        old_miss_hgnc = c(old_miss_hgnc, rep(miss_hgnc[j], length(symbol)))
      } else if (length(prev_symbol) > 0){
        upd_miss_hgnc = c(upd_miss_hgnc, set_hgnc$symbol[prev_symbol])
        old_miss_hgnc = c(old_miss_hgnc, rep(miss_hgnc[j], length(prev_symbol)))
      } else if (length(alias_symbol) > 0){
        upd_miss_hgnc = c(upd_miss_hgnc, set_hgnc$symbol[alias_symbol])
        old_miss_hgnc = c(old_miss_hgnc, rep(miss_hgnc[j], length(alias_symbol)))
      } else {
        upd_miss_hgnc = c(upd_miss_hgnc, NA)
        old_miss_hgnc = c(old_miss_hgnc, miss_hgnc[j])
      }
    }
    
    set1 = genes[-indexes_miss_hgnc,]
    set2 = cbind(old_miss_hgnc, upd_miss_hgnc)
    colnames(set2) = colnames(set1)
    genes = rbind(set1, set2)
  }
  
  #Output
  genes
}
