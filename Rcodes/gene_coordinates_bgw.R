library(readxl)
genomes_bgw = read_excel("genome_version.xlsx", sheet = "Version") 
ncbi_genomes = list.files("./GCF/")
dir.create("./results_gene_coordinates/")

for (i in 1:length(ncbi_genomes)){
  index = gregexpr("_", ncbi_genomes[i])[[1]][2]
  GCF = substr(ncbi_genomes[i], 1, index-1)
  txid = genomes_bgw$txid[genomes_bgw$RefSeq_assembly == GCF]
  table = read.delim(paste0("./GCF/", ncbi_genomes[i]), comment.char = "#", header = F, sep = "\t")
  table = table[table$V3 %in% c("gene", "pseudogene"),]
  table = table[,-c(2,6,8)]
  
  ID = Name = GeneID = c()
  for (j in 1:nrow(table)){
    V9 = unlist(strsplit(table$V9[j], ";"))
    index = grep("^ID=", V9)
    if (length(index) > 0){
      ID[j] = gsub("^ID=", "", V9[index])
    } else {Name[j] = NA}
    
    index = grep("^Name=", V9)
    if (length(index) > 0){
      name = unlist(strsplit(V9[index], "="))[2]
      name = gsub("^.+:", "", name)
      Name[j] = name
    } else {Name[j] = NA}
    
    index = grep("^Dbxref=", V9)
    if (length(index) > 0){
      geneid = gsub("Dbxref=", "", V9[index])
      geneid = unlist(strsplit(geneid, ","))
      geneid = geneid[grep("GeneID:", geneid)]
      geneid = gsub("GeneID:", "", geneid)
      GeneID[j] = geneid
    } else {GeneID[j] = NA}
  }
  table = cbind(ID, Name, GeneID, table)
  
  ncbi_taxon = paste0("NCBITaxon_", txid)
  genome = read_excel("genome_version.xlsx", sheet = ncbi_taxon)
  colnames(genome)[4] = "V1"
  table = merge(table, genome[c(3,4)], by = "V1", all.x = T)
  
  indexes = which(!is.na(table$Name.y))
  table$Name.y[indexes] = paste0("chr", table$Name.y[indexes])
  table$Name.y[-indexes] = table$V1[-indexes]
  
  export_table = table[c(2:4,10,6:8)]
  export_table = export_table[!duplicated(export_table),]
  export_table = cbind(export_table, assembly = genomes_bgw$RefSeq_assembly[genomes_bgw$txid == txid])
  
  indexes = which(is.na(export_table[1]))
  if (length(indexes) > 0){export_table = export_table[-indexes,]}
  
  export_table[is.na(export_table)] = "-"

  write.table(export_table, paste0("./results_gene_coordinates/", ncbi_taxon, ".tsv"), 
              col.names = F, row.names = F, quote = F, sep = "\t")
}
