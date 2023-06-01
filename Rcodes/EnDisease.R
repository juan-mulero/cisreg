#1. Reading and preprocess:
##Linux:
#download.file("http://health.tsinghua.edu.cn/jianglab/endisease/download_files/endisease.tar.gz", 
#              "./EnDisease.tar.gz")
#system("tar xzvf EnDisease.tar.gz")

##Windows:
#download.file("http://health.tsinghua.edu.cn/jianglab/endisease/download_files/endisease.zip", 
#              "./EnDisease.zip")
#unzip("EnDisease.zip")
t = proc.time()
orig_EnDisease = read.delim("endisease_2.0.txt", header = F, sep = "\t", encoding = "UTF-8", stringsAsFactors = F)
EnDisease = orig_EnDisease
colnames(EnDisease) = c("original_ID", "PMID", "disease", "OMIM", "species", "biosample_name", "chr",
                        "start", "end", "strand", "external_gene_name_target_genes", "assembly", "refsnp_ID")

#1.1.Human reference genome:
EnDisease = EnDisease[EnDisease$species == "Homo sapiens",] #Only human enhancers
EnDisease[EnDisease == "/N"] = "-" #representation of null values
indexes = which(EnDisease$assembly != "hg38")
if (length(indexes > 0)) {EnDisease = EnDisease[-indexes,]}

#1.2. Include enhancer coordinates
source("./scripts/bed2coord.R")
coord = bed2coord(EnDisease[7:9])
EnDisease = data.frame(EnDisease[c(1:4,6:9)], coord, EnDisease[11:13])


#2. Enhancer
#enh38
enh38 = EnDisease[c(6:9,11)]
enh38 = dplyr::mutate(enh38, enh_ID = paste0(assembly, "_", chr, "_", start, "_", end))
enh38 = enh38[c(6,1:5)]
enh38 = enh38[!duplicated(enh38),]

#enhancers
enhancers = enh38[-5]
colnames(enhancers)[2:5] = c("orig_chr", "orig_start", "orig_end", "orig_assembly")
enhancers = cbind(enhancers, current_chr = enh38$chr, current_start = enh38$start, current_end = enh38$end, 
                  current_assembly = enh38$assembly, minimum_ratio = "-", score = "-")

#Include enh_ID in the original set
EnDisease = merge(EnDisease, enh38[c(1,5)], by = "coord", all.x = T)
EnDisease = EnDisease[c(13,2:9,1,10:12)]

#Metadata
EnDisease$biosample_name = gsub("\\s+", "", EnDisease$biosample_name)
EnDisease$biosample_name = gsub("/", ",", EnDisease$biosample_name) #The / character is used to separate different values
EnDisease$biosample_name[EnDisease$biosample_name == "BJ,ET"] = "BJ/ET" #But the BJ/ET cell line also uses this character
EnDisease = EnDisease[!duplicated(EnDisease),]
enh_metadata = cbind(EnDisease[1:2], crossref = "-", EnDisease[c(3,6)], enh_method = "-", type = "-",
                     source = "EnDisease")
enh_metadata = enh_metadata[!duplicated(enh_metadata),]
colnames(enh_metadata)[4] = "enh_PMID"

enh_metadata_split = enh_metadata
source("./scripts/SplitIntoAtomicValues.R")
enh_metadata_split = SplitIntoAtomicValues(enh_metadata_split, 5, ",") #split biosamples
enh_metadata_split = enh_metadata_split[!duplicated(enh_metadata_split),]


#3. Target genes
enh2target_genes = cbind(EnDisease[c(1,11,3)], enh2gene_method = "-", EnDisease[6], source = "EnDisease")
colnames(enh2target_genes)[3] = "enh2gene_PMID"
indexes = which(enh2target_genes$external_gene_name_target_genes == "-")
if(length(indexes) > 0){enh2target_genes = enh2target_genes[-indexes,]}
enh2target_genes = SplitIntoAtomicValues(enh2target_genes, 2, "/") #split genes
upd_target_genes = unique(enh2target_genes$external_gene_name_target_genes)
source("./scripts/UpdateGeneSymbol.R")
upd_target_genes = UpdateGeneSymbol(upd_target_genes)
colnames(upd_target_genes) = c("external_gene_name_target_genes", "hgnc_symbol_target_genes")
upd_target_genes[is.na(upd_target_genes)] = "-"
enh2target_genes = merge(enh2target_genes, upd_target_genes, by = "external_gene_name_target_genes", all.x = T)
enh2target_genes = enh2target_genes[!duplicated(enh2target_genes),]
enh2target_genes = enh2target_genes[c(2,1,7,3:6)]

enh2target_genes_split = enh2target_genes
enh2target_genes_split = SplitIntoAtomicValues(enh2target_genes_split, 6, ",") #Split biosamples
enh2target_genes_split = enh2target_genes_split[!duplicated(enh2target_genes_split),]


#4. Diseases
enh2diseases = cbind(EnDisease[c(1,4,5)], DO = "-", mesh = "-", disease_PMID = EnDisease$PMID, disease_method = "-", 
                     source = "EnDisease")
enh2diseases = enh2diseases[!duplicated(enh2diseases),]
indexes = which(enh2diseases$disease == "-")
if (length(indexes) > 0){enh2diseases = enh2diseases[-indexes,]}
indexes = which(grepl("^[[:alnum:][:blank:][:punct:]]+$", enh2diseases$disease) == F)
#unique(enh2diseases$disease[indexes]) #"Bazex-Dupr??????-Christol syndrome"
if (length(indexes) > 0){enh2diseases$disease[indexes] = "-"}
enh2diseases$disease = trimws(enh2diseases$disease)


#5. Mutations
enh2mutations = cbind(EnDisease[c(1,13)], mutation_PMID = EnDisease$PMID, mutation_method = "-", source = "EnDisease")
indexes = which(enh2mutations$refsnp_ID == "-")
if(length(indexes) > 0){enh2mutations = enh2mutations[-indexes,]}
indexes = grep("^rs", enh2mutations$refsnp_ID)
#enh2mutations$refsnp_ID #"rs2435357"  "rs3858145"  "rs2238126"  "rs12913832" "4q32 A>C"   "rs2736100" 
enh2mutations = enh2mutations[indexes,]


#6. Biosamples
biosample_name = c(enh_metadata_split$biosample_name, enh2target_genes_split$biosample_name)
biosample_name = unique(biosample_name)
biosample_name = data.frame(biosample_name)
index = which(biosample_name$biosample_name == "-")
if(length(index) > 0){biosample_name = biosample_name[-index,]}
biosample_name = data.frame(biosample_name) #table for mapping terms in ontologies


#7. Merge data in an unique file
rm(coord, index, indexes)
#At the moment we have different subsets of data according to different topics. Now we join all these subsets in a unified format
merge_EnDisease = merge(enhancers, enh_metadata_split, by = "enh_ID", all = T)
merge_EnDisease = merge(merge_EnDisease, enh2target_genes_split[1:5], by = "enh_ID", all = T)
merge_EnDisease = merge_EnDisease[!duplicated(merge_EnDisease),]
merge_EnDisease = cbind(merge_EnDisease, external_gene_name_TFs = "-", hgnc_symbol_TFs = "-", TFs2enh_PMID = "-",
                        TFs2enh_method = "-")
merge_EnDisease = merge(merge_EnDisease, enh2diseases[1:7], by = "enh_ID", all = T)
merge_EnDisease = merge_EnDisease[!duplicated(merge_EnDisease),]
merge_EnDisease = merge(merge_EnDisease, enh2mutations[1:4], by = "enh_ID", all = T)
merge_EnDisease = merge_EnDisease[!duplicated(merge_EnDisease),]
merge_EnDisease[is.na(merge_EnDisease)] = "-"
merge_EnDisease$original_ID = "-" #because the identifier is not specific to the enhancer, but to the enhancer-disease link.

#We remove the labels that we are not going to use
filt_merge_EnDisease = merge_EnDisease[-c(19,23,27)] #external_gene_name and disease
filt_merge_EnDisease = filt_merge_EnDisease[!duplicated(filt_merge_EnDisease),]
indexes = which(filt_merge_EnDisease$hgnc_symbol_target_genes == "-")
if (length(indexes) > 0) {filt_merge_EnDisease$enh2gene_PMID[indexes] = filt_merge_EnDisease$enh2gene_method[indexes] = "-"}
indexes = which(filt_merge_EnDisease$OMIM == "-")
if (length(indexes) > 0) {filt_merge_EnDisease$disease_PMID[indexes] = filt_merge_EnDisease$disease_method[indexes] = "-"}
indexes = which(filt_merge_EnDisease$refsnp_ID == "-")
if (length(indexes) > 0) {filt_merge_EnDisease$mutation_PMID[indexes] = filt_merge_EnDisease$mutation_method[indexes] = "-"}

#Merge diseases: We represent the disease identifiers in a single variable
disease = c()
for (i in 1:nrow(filt_merge_EnDisease)){
  OMIM = filt_merge_EnDisease$OMIM[i]
  if (OMIM != "-"){OMIM = paste0("OMIM:", OMIM)}
  mesh = filt_merge_EnDisease$mesh[i]
  if (mesh != "-"){mesh = paste0("mesh:", mesh)}
  DO = filt_merge_EnDisease$DO[i]
  list = c(OMIM, DO, mesh)
  list = unique(list)
  if (length(list) > 1 & "-" %in% list){
    index  = which(list == "-")
    list = list[-index]
  }
  disease[i] = paste(list, collapse = ",")
}

filt_merge_EnDisease = cbind(filt_merge_EnDisease[1:24], disease, filt_merge_EnDisease[28:32])
filt_merge_EnDisease = SplitIntoAtomicValues(filt_merge_EnDisease, 25, ",") #Split diseases
filt_merge_EnDisease = filt_merge_EnDisease[!duplicated(filt_merge_EnDisease),]
filt_merge_EnDisease$original_ID = "-"
indexes = which(filt_merge_EnDisease$disease == "-")
if (length(indexes) > 0) {filt_merge_EnDisease$disease_PMID[indexes] = filt_merge_EnDisease$disease_method[indexes] = "-"}


#8. Save
dir.create("./EnDisease_results")
rm(disease, DO, i, index, indexes, list, mesh, OMIM)
t = proc.time() - t; t
save.image("./EnDisease_results/Data_EnDisease.RData")
write.table(filt_merge_EnDisease, "./EnDisease_results/EnDisease.tsv", quote = F, sep = "\t", col.names = T, row.names = F, 
            fileEncoding = "utf8")
