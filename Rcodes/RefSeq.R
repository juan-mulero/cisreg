#Source: https://www.ncbi.nlm.nih.gov/genome/guide/human/
#Execution: 30/12/2021
#annotation-date 11/19/2021

#1. Reading files:
library(data.table)
download.file("https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz", 
              "GRCh38_latest_genomic.gff.gz")
#system("tar xzvf GRCh38_latest_genomic.gff.gz")
hg38 = fread("GRCh38_latest_genomic.gff", header = F, sep = "\t", skip = 9, fill = T)
#sort(unique(hg38$V3))
orig_enh38 = hg38[hg38$V3 == "enhancer",]
#dim(orig_enh38) #5643    9
#table(enhancers$V2)
#Curated Genomic        RefSeqFE 
#8            5635

##Chromosomes are identified as LOCUS (example: NC_000001). We modify this format to chrN
enhancers = orig_enh38
indexes = grep("NC_",enhancers$V1)
#length(indexes) #5326
enhancers = enhancers[indexes,]
enhancers$V1 = as.integer(substr(enhancers$V1,8,9))
enhancers$V1[enhancers$V1 == 23] = "X"
enhancers$V1[enhancers$V1 == 24] = "Y"
enhancers$V1 = paste0("chr", enhancers$V1)


#2. Data frame with the attributes
##First we identify all attribute keys
key = as.character(enhancers[1,9])
key = unlist(strsplit(key, split = ";"))
key = as.character(key[1])
key = unlist(strsplit(key, split = "="))
key = as.character(key[1])

keys = c(key)
keys_attributes = data.frame(keys)
  
for (i in 1:nrow(enhancers)){
  attributes = as.character(enhancers[i,9])
  attributes = strsplit(attributes, split = ";")
  attributes = unlist(attributes)
  keys = c()
  
  for (j in 1:length(attributes)){
    feature = as.character(attributes[j])
    feature = strsplit(feature, split = "=")
    feature = unlist(feature)
    keys[j] = as.character(feature[1])
  }
  keys_attributes = merge(data.frame(keys), keys_attributes, all = T)
}
  
##Then, we fill in the table:
df_attributes = matrix(NA, nrow = nrow(enhancers), ncol = nrow(keys_attributes))
colnames(df_attributes) = keys_attributes$keys

for (i in 1:nrow(enhancers)){
  attributes = as.character(enhancers[i,9])
  attributes = strsplit(attributes, split = ";")
  attributes = unlist(attributes)
  
  for (j in 1:length(attributes)){
    feature = as.character(attributes[j])
    feature = strsplit(feature, split = "=")
    feature = unlist(feature)
    key = as.character(feature[1])
    value = as.character(feature[2])
    
    df_attributes[i, key] = value
  }
}
#dim(df_attributes) #5326    9
df_attributes = data.frame(df_attributes)
rm(keys_attributes, attributes, feature, i, indexes, j, key, keys, value)


#3. Enhancers
enhancers = data.frame(enhancers)
enhancers = cbind(enhancers[c(1:5,7)], df_attributes[c(1:3,6:7)])

##Processing of the column Dbxref
source("./scripts/SplitIntoAtomicValues.R")
enhancers = SplitIntoAtomicValues(enhancers, 7, ",")
keys = unique(unlist(lapply(enhancers$Dbxref, function(ref){unlist(strsplit(ref, ":"))[1]})))
#keys #"GeneID" "VISTA"  "HGNC"   "MIM"    

###HGNC to HGNC crossref of HUGO
enhancers$Dbxref = gsub("HGNC:HGNC:", "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/", enhancers$Dbxref)

###VISTA to VISTA crossref
indexes = grep("VISTA", enhancers$Dbxref)
enhancers$Dbxref = gsub("VISTA:hs", "", enhancers$Dbxref)
enhancers$Dbxref[indexes] = paste0("https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?form=presentation&show=1912&experiment_id=",
                                   enhancers$Dbxref[indexes], "&organism_id=1")
###GeneID to NCBI gene crossref
indexes = grep("GeneID:", enhancers$Dbxref)
enhancers$Dbxref = gsub("GeneID:", "", enhancers$Dbxref)
enhancers$Dbxref[indexes] = paste0("https://www.ncbi.nlm.nih.gov/gene/", enhancers$Dbxref[indexes])

###MIM to OMIM crossref
indexes = grep("MIM:", enhancers$Dbxref)
enhancers$Dbxref = gsub("MIM:", "", enhancers$Dbxref)
enhancers$Dbxref[indexes] = paste0("https://www.omim.org/entry/", enhancers$Dbxref[indexes])


##Processing of the column experiment
enhancers = SplitIntoAtomicValues(enhancers, 8, ",")
enh_method = enh_method_ont = enh_PMID = c()
for (i in 1:nrow(enhancers)){
  experiment = enhancers$experiment[i]
  if(!is.na(experiment)){
    split = unlist(strsplit(experiment, "\\["))
    split = trimws(split)
    split = gsub("\\]", "", split)
    enh_method[i] = split[1]
    enh_method_ont[i] = split[2]
    enh_PMID[i] = split[3]
  }
}
enhancers = cbind(enhancers[1:7], enh_method, enh_method_ont, enh_PMID, enhancers[9:11])
enhancers = SplitIntoAtomicValues(enhancers, 10, " ")
enhancers$enh_PMID = gsub("%2C", "", enhancers$enh_PMID)
enhancers$enh_PMID = gsub("PMID:", "", enhancers$enh_PMID)

enhancers$enh_method = unlist(lapply(enhancers$enh_method, function(i){unlist(strsplit(i, ":"))[2]}))


##Processing of the column ID
enhancers$ID = gsub("id-", "", enhancers$ID)


##Final table
source("./scripts/bed2coord.R")
coord = bed2coord(enhancers[c(1,4,5)])
refseq_enh = cbind(enh_ID = paste0("hg38_", gsub("[:-]", "_", coord)), orig_chr = enhancers$V1, orig_start = enhancers$V4,
                   orig_end = enhancers$V5, orig_assembly = "hg38", current_chr = enhancers$V1, current_start = enhancers$V4,
                   current_end = enhancers$V5, current_assembly = "hg38", minimum_ratio = NA, score = NA, original_ID = enhancers$ID,
                   crossref = enhancers$Dbxref, enhancers[10], biosample_name = NA, enhancers[8], type = enhancers$V3, 
                   source = "RefSeq", hgnc_symbol_target_genes = NA, enh2gene_PMID = NA, enh2gene_method = NA, hgnc_symbol_TFs = NA,
                   TFs2enh_PMID = NA, TFs2enh_method = NA, disease = NA, disease_PMID = NA, disease_method = NA, refsnp_ID = NA,
                   mutation_PMID = NA, mutation_method = NA)
refseq_enh[is.na(refseq_enh)] = "-"

#4. Save files
rm(coord, enh_method, enh_method_ont, enh_PMID, experiment, i, indexes, keys, split)
dir.create("./RefSeq_results")
save.image("./RefSeq_results/Data_RefSeq.RData")
write.table(refseq_enh, "./RefSeq_results/RefSeq.tsv", quote = F, sep = "\t", col.names = T, row.names = F, 
            fileEncoding = "utf8")
