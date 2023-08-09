#Code to extract from HGNC the current gene symbol corresponding to one alias or to a previous gene symbol. 
#For this purpose, REST queries are used.
#Input: vectorGenes --> vector with the genes.
#Output: vectorGenes --> vector with the updated genes.
queryHGNC = function(vectorGenes){
  library(xml2)
  library(rvest)
  
  ExtractSymbol = function(url){
    read = read_xml(url)
    node = html_nodes(read, xpath = '//str[contains(@name,"symbol")]')
    node = xml_text(node)
    node
  }

  for (i in 1:length(vectorGenes)){
    url = paste0("http://rest.genenames.org/fetch/alias_symbol/", vectorGenes[i])
    symbol = ExtractSymbol(url)
    if (length(symbol) == 1) {gene = symbol}
    else if (length(symbol) > 1){gene = NA} 
    else {
      url = paste0("http://rest.genenames.org/fetch/prev_symbol/", vectorGenes[i])
      symbol = ExtractSymbol(url)
      if (length(symbol) == 1) {gene = symbol}
      else if (length(symbol) > 1) {gene = NA} 
      else { 
        url = paste0("http://rest.genenames.org/fetch/symbol/", vectorGenes[i])
        symbol = ExtractSymbol(url)
        if (length(symbol) == 1) {gene = symbol}
        else {gene = NA}
      }
    }
    vectorGenes[i] = gene
  }
  vectorGenes
}