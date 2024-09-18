#Function to obtain the OMIM identifiers linked to Disease Ontology through the DO ontology.
#Input: vector with DOID
#Output: vector with OMIM
DO2OMIM = function(DOID){
  DiseaseOntology = read.delim("https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/HumanDO.obo", 
                               header = F, sep = "\n")
  term = which(DiseaseOntology$V1 == "[Term]")
  output = c()
  for (i in 1:length(DOID)){
    if (!is.na(DOID[i]) & DOID[i] != "-"){
      if (grepl("^DOID:", DOID[i])){
        index = which(DiseaseOntology$V1 == paste0("id: ", DOID[i]))
      } 
      else {index = which(DiseaseOntology$V1 == paste0("id: DOID:", DOID[i]))}
      
      if (length(index) > 0){
        data = DiseaseOntology$V1[(index-1):(term[which(term > index)[1]]-1)]
        indexes = grep("^xref: O?MIM:", data)
        if (length(indexes) > 0){
          values = c()
          for (j in indexes){
            value = unlist(strsplit(data[j], "O?MIM:"))[2]
            values = c(values, value)
          }
          output[i] = paste(values, collapse = ",")
        } else {output[i] = NA}
      }
    } else {output[i] = NA}
  }
  output
}
