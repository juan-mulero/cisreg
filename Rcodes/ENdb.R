#1. Reading and preprocess:
orig_ENdb = read.delim("http://www.licpathway.net/ENdb/file/download/ENdb_enhancer.txt", 
                           header = T, sep = "\t", encoding = "UTF-8")
orig_ENdb = data.frame(orig_ENdb)

#1.1.Human reference genome:
hg_ENdb = orig_ENdb
hg_ENdb[hg_ENdb == "--" | hg_ENdb == "" | hg_ENdb == " " | is.na(hg_ENdb)] = "-" #Different ways of representing null values
#unique(hg_ENdb$Reference_genome) #"hg19" "mm10"
hg_ENdb = hg_ENdb[hg_ENdb$Reference_genome == "hg19",] #human enhancers

#1.2.Check of chr
#unique(sort(hg_ENdb$Chromosome))
#[1] "ch1"   "ch3"   "chr1"  "chr10" "chr11" "Chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr2"  "chr20"
#[17] "chr21" "chr22" "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9"  "chr9 " "chrx"  "chrX" 
hg_ENdb$Chromosome = stringr::str_to_lower(hg_ENdb$Chromosome)
hg_ENdb$Chromosome = gsub("x", "X", hg_ENdb$Chromosome)
hg_ENdb$Chromosome = gsub("ch", "chr", hg_ENdb$Chromosome)
hg_ENdb$Chromosome = gsub("chrr", "chr", hg_ENdb$Chromosome)
hg_ENdb$Chromosome = gsub("\\s", "", hg_ENdb$Chromosome)

#1.3.Include enhancer coordinates
source("./scripts/bed2coord.R")
coord = bed2coord(hg_ENdb[5:7])
hg_ENdb = cbind(hg_ENdb[c(1,2,5:7)], coord, hg_ENdb[c(9,11,13:16,19:21,25,27,29,31)])
rm(coord)


#2. Enhancers:
#2.1.hg19
enh19 = hg_ENdb[3:6]
colnames(enh19) = c("orig_chr", "orig_start", "orig_end", "coord")
enh19 = enh19[!duplicated(enh19),]
enh19 = cbind(enh19, orig_assembly = "hg19")
enh19 = dplyr::mutate(enh19, enh_ID = paste0(orig_assembly, "_", orig_chr, "_", orig_start, "_", orig_end))
enh19 = enh19[c(6,1:5)]

#2.2.hg38
source("./scripts/hgtohg_onestep.R")
hg19tohg38 = hgtohg_onestep(enh19[2:6], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
hg19tohg38_tables = hgtohg2tables("hg19", "hg38", enh19, hg19tohg38)


#2.3.Enhancer metadata
#ID:
hg_ENdb = merge(hg_ENdb, enh19[c(1,5)], by = "coord", all.x = T)
hg_ENdb = hg_ENdb[c(20,2:6,1,7:19)]

#Crossref:
crossref = unlist(lapply(hg_ENdb$Enhancer_id, function(ID){
  paste0("http://www.licpathway.net/ENdb/search/Detail.php?Species=Human&Enhancer_id=", ID)}))
hg_ENdb = cbind(hg_ENdb[1:2], crossref, hg_ENdb[3:20])

#Type:
#unique(hg_ENdb$Enhancer_type) #[1] "Enhancer"       "Super-Enhancer" "enhancer"       "Enhancer "  
hg_ENdb$Enhancer_type = stringr::str_to_lower(hg_ENdb$Enhancer_type)
hg_ENdb$Enhancer_type = gsub("\\s", "", hg_ENdb$Enhancer_type)

#Source:
hg_ENdb = cbind(hg_ENdb, source = "ENdb")

#Set of metadata:
enh_metadata = hg_ENdb[c(1:4,9:11,22)]
colnames(enh_metadata)[2:7] = c("original_ID", "crossref", "enh_PMID", "biosample_name", "enh_method", "type")
enh_metadata_split = enh_metadata
source("./scripts/SplitIntoAtomicValues.R")
enh_metadata_split = SplitIntoAtomicValues(enh_metadata_split, 5, ",")
enh_metadata_split = SplitIntoAtomicValues(enh_metadata_split, 6, ",")
enh_metadata_split = enh_metadata_split[!duplicated(enh_metadata_split),]
enh_metadata_split$biosample_name = trimws(enh_metadata_split$biosample_name)
enh_metadata_split$enh_method = trimws(enh_metadata_split$enh_method)
enh_metadata_split = enh_metadata_split[!duplicated(enh_metadata_split),]
indexes = which(grepl("^[[:alnum:][:blank:][:punct:]]+$", enh_metadata_split$biosample_name) == F) #to remove cell lines with strange characters
#unique(enh_metadata_split$biosample_name[indexes])
#[1] "Immortalized Murine Gonadotrope L????????????T2 Cell" "231/ER????????????Wt Cells"                          
#[3] "231/ER????????????LQ Cells"                           "????????????-Cell"                                   
#[5] "????????????-cell"                                    "????????????TC-3"                                    
#[7] "????????????TC-6" 
if (length(indexes) > 0){enh_metadata_split$biosample_name[indexes] = "-"}


#3. Target genes
enh2target_genes = hg_ENdb[c(1,12,4,13,14,9,22)]

#Methods:
enh2gene_method = c()
for (i in 1:nrow(enh2target_genes)){
  strong_exp = enh2target_genes$target_gene_strong_experiment[i]
  weak_exp = enh2target_genes$target_gene_weak_experiment[i]
  if (strong_exp != "-" & weak_exp != "-")
  {enh2gene_method[i] = paste0(strong_exp, ",", weak_exp)}
  else if (strong_exp != "-" & weak_exp == "-"){
    enh2gene_method[i] = strong_exp}
  else {enh2gene_method[i] = weak_exp}
}
enh2target_genes = cbind(enh2target_genes[1:3], enh2gene_method, enh2target_genes[6:7])

#Gene symbol update:
enh2target_genes = SplitIntoAtomicValues(enh2target_genes, 2, ",") #Split genes
colnames(enh2target_genes)[c(2,3,5)] = c("external_gene_name", "enh2gene_PMID", "biosample_name")
enh2target_genes$external_gene_name = trimws(enh2target_genes$external_gene_name)
source("./scripts/UpdateGeneSymbol.R")
upd_target_genes = UpdateGeneSymbol(enh2target_genes$external_gene_name)
enh2target_genes = merge(enh2target_genes, upd_target_genes, by = "external_gene_name")
enh2target_genes = enh2target_genes[c(2,1,7,3:6)]
enh2target_genes = enh2target_genes[!duplicated(enh2target_genes),]

#Set of metadata:
enh2target_genes$enh2gene_method = trimws(enh2target_genes$enh2gene_method)
try({enh2target_genes$biosample_name = trimws(enh2target_genes$biosample_name)}, silent = T) #due to strange characters
enh2target_genes_split = enh2target_genes
enh2target_genes_split = SplitIntoAtomicValues(enh2target_genes_split, 5, ",") #Split the methods
enh2target_genes_split$enh2gene_method = trimws(enh2target_genes_split$enh2gene_method)
enh2target_genes_split = SplitIntoAtomicValues(enh2target_genes_split, 6, ",") #Split the biosamples
enh2target_genes_split$biosample_name = trimws(enh2target_genes_split$biosample_name)
enh2target_genes_split = enh2target_genes_split[!duplicated(enh2target_genes_split),]
indexes = which(grepl("^[[:alnum:][:blank:][:punct:]]+$", enh2target_genes_split$biosample_name) == F)
if (length(indexes) > 0){enh2target_genes_split$biosample_name[indexes] = "-"}


#4. Diseases
enh2diseases = hg_ENdb[c(1,15,16,17,4,22)]
enh2diseases$Disease = trimws(enh2diseases$Disease)
indexes = which(enh2diseases$Disease == "-")
if(length(indexes) > 0){enh2diseases = enh2diseases[-indexes,]}
enh2diseases = enh2diseases[!duplicated(enh2diseases),]
colnames(enh2diseases)[c(2,5)] = c("disease", "disease_PMID")
enh2diseases = cbind(enh2diseases[1:2], OMIM = "-", enh2diseases[3:5], disease_method = "-", enh2diseases[6])

#unique(enh2diseases$disease)
#When we check the data, we find incorrect annotations that we remove from the dataset
indexes = which(enh2diseases$disease %in% c("National Institute of Allergy and Infectious Diseases", 
                                            "Mycobacterium Tuberculosis", "Skin Barrier"))
if (length(indexes > 0)) {enh2diseases = enh2diseases[-indexes,]}

#We remove multiple records (with commas) because there is no direct association between names and identifiers,
#but we also create a file with this data adapted.
other_enh2diseases = enh2diseases
indexes = grep(",", enh2diseases$disease)
if (length(indexes) > 0) {
  enh2diseases = enh2diseases[-indexes,]
  
  count = stringr::str_count(other_enh2diseases$disease[indexes], ",") + 1 
  disease_name = paste(other_enh2diseases$disease[indexes], collapse = ",")
  disease_name = unlist(strsplit(disease_name, ","))
  data = data.frame(cbind(enh_ID = rep(other_enh2diseases$enh_ID[indexes], count), disease = disease_name, OMIM = "-",  
               DO = c("-", "DOID:417", "-", "DOID:9256", "DOID:1612", "DOID:3393", "DOID:6364", "DOID:9348", "-", "DOID:10763"),
               mesh = c("D006084", "D001327", "D000022", "D015179", "D001943", "D003324", "D008881", "-", "D005352", "D006973"), 
               disease_PMID = rep(other_enh2diseases$disease_PMID[indexes], count), disease_method = ",", source = "ENdb"))
  other_enh2diseases = other_enh2diseases[-indexes,]
  other_enh2diseases = rbind(other_enh2diseases, data)
}

#Mapping with OMIM
diseases = other_enh2diseases[2:5]
diseases = diseases[!duplicated(diseases),]
source("./scripts/DO2OMIM.R")
diseases$OMIM = DO2OMIM(diseases$DO)
diseases = SplitIntoAtomicValues(diseases, 2, ",")
diseases[is.na(diseases)] = "-"
diseases = diseases[!duplicated(diseases),]

other_enh2diseases = merge(other_enh2diseases[-3], diseases[1:2], by = "disease", all = T)
other_enh2diseases = other_enh2diseases[c(2,1,8,3:7)]


#5. Transcription factors
TFs2enh = hg_ENdb[c(1,18,4,19,9,22)]
indexes = which(TFs2enh$TF_name == "-")
if (length(indexes) > 0){TFs2enh = TFs2enh[-indexes,]}
TFs2enh$TF_name = trimws(TFs2enh$TF_name)
TFs2enh = SplitIntoAtomicValues(TFs2enh, 2, ",") #Split TFs
colnames(TFs2enh)[2:5] = c("external_gene_name", "TFs2enh_PMID", "TFs2enh_method", "biosample_name")

#Gene symbol update
upd_TFs = UpdateGeneSymbol(TFs2enh$external_gene_name)
TFs2enh = merge(TFs2enh, upd_TFs, by = "external_gene_name")
TFs2enh = TFs2enh[c(2,1,7,3:6)]
TFs2enh = TFs2enh[!duplicated(TFs2enh),]

#Set of metadata
TFs2enh$TFs2enh_method = trimws(TFs2enh$TFs2enh_method)
try({TFs2enh$biosample_name = trimws(TFs2enh$biosample_name)}, silent = T)
TFs2enh_split = TFs2enh
TFs2enh_split = SplitIntoAtomicValues(TFs2enh_split, 5, ",") #Split methods
TFs2enh_split$TFs2enh_method = trimws(TFs2enh_split$TFs2enh_method)
TFs2enh_split = SplitIntoAtomicValues(TFs2enh_split, 6, ",") #Biosamples
TFs2enh_split$biosample_name = trimws(TFs2enh_split$biosample_name)
TFs2enh_split = TFs2enh_split[!duplicated(TFs2enh_split),]
indexes = which(grepl("^[[:alnum:][:blank:][:punct:]]+$", TFs2enh_split$biosample_name) == F)
if (length(indexes) > 0){TFs2enh_split$biosample_name[indexes] = "-"} 


#6. Mutations
enh2mutations = hg_ENdb[c(1,20,4,21,22)]
indexes = which(enh2mutations$SNP_id == "-")
if (length(indexes > 0)){enh2mutations = enh2mutations[-indexes,]}
enh2mutations = enh2mutations[!duplicated(enh2mutations),]
enh2mutations = SplitIntoAtomicValues(enh2mutations, 2, ",") #Split mutations
enh2mutations$SNP_id = trimws(enh2mutations$SNP_id)
colnames(enh2mutations)[2:4] = c("refsnp_ID", "mutation_PMID", "mutation_method")
enh2mutations = enh2mutations[!duplicated(enh2mutations),]

enh2mutations_split = enh2mutations
enh2mutations_split = SplitIntoAtomicValues(enh2mutations_split, 4, ",") #Spit methods
enh2mutations_split$mutation_method = trimws(enh2mutations_split$mutation_method)
enh2mutations_split = enh2mutations_split[!duplicated(enh2mutations_split),]


#7. Methods
methods = c(enh_metadata_split$enh_method, enh2target_genes_split$enh2gene_method, TFs2enh_split$TFs2enh_method, 
           enh2mutations_split$mutation_method)
methods = unique(methods)
methods = data.frame(methods)
index = which(methods$methods == "-")
if(length(index) > 0){methods = methods[-index,]}
methods = data.frame(methods) #table for mapping terms in ontologies


#8. Biosamples
biosample_names = c(enh_metadata_split$biosample_name, enh2target_genes_split$biosample_name, TFs2enh_split$biosample_name)
biosample_names = unique(biosample_names)
biosample_names = data.frame(biosample_names)
index = which(biosample_names$biosample_names == "-")
if(length(index) > 0){biosample_names = biosample_names[-index,]}
biosample_names = data.frame(biosample_names) #table for mapping terms in ontologies


#9. Merge data
rm(count, crossref, disease_name, enh2gene_method, i, index, indexes, strong_exp, weak_exp)
#At the moment we have different subsets of data according to different topics. Now we join all these subsets in an unified format
merge_ENdb = merge(hg19tohg38_tables$enhancers, enh_metadata_split, by = "enh_ID", all = T)
merge_ENdb = merge(merge_ENdb, enh2target_genes_split[1:5], by = "enh_ID", all = T)
merge_ENdb = merge_ENdb[!duplicated(merge_ENdb),]
colnames(merge_ENdb)[19:20] = c("external_gene_name_target_genes", "hgnc_symbol_target_genes")
merge_ENdb = merge(merge_ENdb, TFs2enh_split[1:5], by = "enh_ID", all = T)
merge_ENdb = merge_ENdb[!duplicated(merge_ENdb),]
colnames(merge_ENdb)[23:24] = c("external_gene_name_TFs", "hgnc_symbol_TFs")
merge_ENdb = merge(merge_ENdb, enh2diseases[1:7], by = "enh_ID", all = T)
merge_ENdb = merge_ENdb[!duplicated(merge_ENdb),]
merge_ENdb = merge(merge_ENdb, enh2mutations_split[1:4], by = "enh_ID", all = T)
merge_ENdb = merge_ENdb[!duplicated(merge_ENdb),]
merge_ENdb[is.na(merge_ENdb)] = "-"

#We remove the labels that we are not going to use
filt_merge_ENdb = merge_ENdb[-c(19,23,27)] #external_gene_name and disease
filt_merge_ENdb = filt_merge_ENdb[!duplicated(filt_merge_ENdb),]
indexes = which(filt_merge_ENdb$hgnc_symbol_target_genes == "-")
if (length(indexes) > 0) {filt_merge_ENdb$enh2gene_PMID[indexes] = filt_merge_ENdb$enh2gene_method[indexes] = "-"}
indexes = which(filt_merge_ENdb$hgnc_symbol_TFs == "-")
if (length(indexes) > 0) {filt_merge_ENdb$TFs2enh_PMID[indexes] = filt_merge_ENdb$TFs2enh_method[indexes] = "-"}
indexes = which(filt_merge_ENdb$refsnp_ID == "-")
if (length(indexes) > 0) {filt_merge_ENdb$mutation_PMID[indexes] = filt_merge_ENdb$mutation_method[indexes] = "-"}

#We represent the disease identifiers in a single variable
disease = c()
for (i in 1:nrow(filt_merge_ENdb)){
  OMIM = filt_merge_ENdb$OMIM[i]
  if (OMIM != "-"){OMIM = paste0("OMIM:", OMIM)}
  mesh = filt_merge_ENdb$mesh[i]
  if (mesh != "-"){mesh = paste0("mesh:", mesh)}
  DO = filt_merge_ENdb$DO[i]
  list = c(OMIM, DO, mesh)
  list = unique(list)
  if (length(list) > 1 & "-" %in% list){
    index  = which(list == "-")
    list = list[-index]
  }
  disease[i] = paste(list, collapse = ",")
}

filt_merge_ENdb = cbind(filt_merge_ENdb[1:24], disease, filt_merge_ENdb[28:32])
filt_merge_ENdb = SplitIntoAtomicValues(filt_merge_ENdb, 25, ",") #Split the diseases
filt_merge_ENdb = filt_merge_ENdb[!duplicated(filt_merge_ENdb),]
indexes = which(filt_merge_ENdb$disease == "-")
if (length(indexes) > 0) {filt_merge_ENdb$disease_PMID[indexes] = filt_merge_ENdb$disease_method[indexes] = "-"}


#10. Save
dir.create("./ENdb_results")
rm(disease, DO, i, index, indexes, list, mesh, OMIM)
save.image("./ENdb_results/Data_ENdb.RData")
write.table(filt_merge_ENdb, "./ENdb_results/ENdb.tsv", quote = F, sep = "\t", col.names = T, row.names = F, 
            fileEncoding = "utf8")

indexes = grep(" and ", filt_merge_ENdb$biosample_name)
set1 = filt_merge_ENdb[indexes,]
set2 = filt_merge_ENdb[-indexes,]
set1 = SplitIntoAtomicValues(set1, 15, " and ")
set1$biosample_name = trimws(set1$biosample_name)
filt_merge_ENdb = rbind(set1, set2)
filt_merge_ENdb = filt_merge_ENdb[order(filt_merge_ENdb$current_chr, filt_merge_ENdb$current_start, filt_merge_ENdb$current_end),]
save.image("./ENdb_results/Data_ENdb.RData")
write.table(filt_merge_ENdb, "./ENdb_results/ENdb.tsv", quote = F, sep = "\t", col.names = T, row.names = F, 
            fileEncoding = "utf8")
