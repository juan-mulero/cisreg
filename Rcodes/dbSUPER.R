#1. Reading:
library(data.table)
orig_dbSUPER = data.frame(fread("https://asntech.org/dbsuper/data/dbSUPER_SuperEnhancers_hg19.tsv"))
dbSUPER = orig_dbSUPER
crossref = unlist(lapply(dbSUPER$se_id, function(ID){
  paste0("https://asntech.org/dbsuper/details.php?se_id=", ID)}))
dbSUPER = cbind(dbSUPER, crossref)
dbSUPER[dbSUPER == ""] = NA #Because there are records without target genes

#Include enhancer coordinates
source("./scripts/bed2coord.R")
coord = bed2coord(dbSUPER[1:3])
dbSUPER = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", coord)), dbSUPER[1:3], coord, dbSUPER[4:8])
rm(coord, crossref)

#Include methods and PMID
##The dataset available for download does not contain this information, but table S1 does. 
##Therefore, we can use this information to enrich the dataset (Table in the ref PMID: 26438538)
cell_name = c("HBL1", "Ly1", "DHL6", "Toledo", "Ly4", "Ly3", "Tonsil", "NCI-H82", "GLC16", "NCI-H69",
              "HSMM", "NHEK", "Pancreatic islets", "H9", "TC32", "TC71")
PMID = c(rep(24332044, 7), rep(25490451, 3), rep(26438538,4), rep(26337082, 2))
data = data.frame(cell_name, PMID)
data = merge(data, dbSUPER[c(6,8)], by = "cell_name", all = T)
data$number = unlist(lapply(data$se_id, function(num){unlist(strsplit(num, "_"))[2]}))
data$number = as.numeric(data$number)
data$method = "H3K27ac" #initialisation with the value most common
indexes = which(is.na(data$PMID))
data$PMID[indexes] = 24119843

##Some cell lines have enhancers from more than one source and with methods different to H3K27ac
data$method[data$cell_name == "Ly1" & data$number >= 65951] = "BRD4"
data$method[data$cell_name == "u87" & data$number >= 67446] = "MED1"
data$PMID[data$cell_name == "u87" & data$number >= 67446] = 23582323
data$PMID[data$cell_name == "H2171" & data$number >= 66864] = 23582323
data$method[data$cell_name == "H2171" & data$number >= 66864] = "MED1"
data$PMID[data$cell_name == "MM1S" & data$number >= 67138] = 23582323
data$method[data$cell_name == "MM1S" & data$number >= 67138] = "MED1"
data$PMID[data$cell_name == "Jurkat" & data$number >= 66235] = 25043025

data = data[-4] #Because we only need the number to construct the IDs
  

#2. Enhancers:
enhancer_metadata = merge(dbSUPER[-c(7,9)], data[2:4], by = "se_id")
colnames(enhancer_metadata)[c(1,3:5,7,9,10)] = c("original_ID", "orig_chr", "orig_start", "orig_end", "biosample_name", 
                                         "enh_PMID", "enh_method")
enhancer_metadata = cbind(enhancer_metadata[c(2:6)], orig_assembly = "hg19", enhancer_metadata[c(1,8,9,7,10)], 
                           type = "super-enhancer", source = "dbSUPER")

#Remap hg38
source("./scripts/hgtohg_onestep.R")
enh19 = enhancer_metadata[1:6]
enh19 = enh19[!duplicated(enh19),]

hg19tohg38 = hgtohg_onestep(enh19[2:5], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
hg19tohg38_tables = hgtohg2tables("hg19", "hg38", enh19, hg19tohg38)
enhancer_metadata = merge(hg19tohg38_tables$enhancers, enhancer_metadata[c(1,7:13)], by = "enh_ID", all = T)
  

#3. Target genes
target_genes = dbSUPER$gene_symbol
target_genes = unique(target_genes)
index = which(is.na(target_genes))
if(length(index) > 0){target_genes = target_genes[-index]}
source("./scripts/UpdateGeneSymbol.R")
target_genes = UpdateGeneSymbol(target_genes)
target_genes = target_genes[!duplicated(target_genes),]

target_genes_metadatada = dbSUPER[c(1,6,7)]
colnames(target_genes_metadatada)[2:3] = c("original_ID", "external_gene_name")
target_genes_metadatada = merge(target_genes_metadatada, target_genes, by = "external_gene_name", all.x = T)
target_genes_metadatada = merge(target_genes_metadatada, enhancer_metadata[c(1,12,14)], by = c("enh_ID", "original_ID"), all.x = T)
target_genes_metadatada = target_genes_metadatada[!duplicated(target_genes_metadatada),]
target_genes_metadatada = cbind(target_genes_metadatada, enh2gene_method = "closest")
colnames(target_genes_metadatada)[5] = "enh2gene_PMID"

merge_dbSUPER = merge(enhancer_metadata, target_genes_metadatada, by = c("enh_ID", "original_ID"), all = T)
merge_dbSUPER = merge_dbSUPER[!duplicated(merge_dbSUPER),]
colnames(merge_dbSUPER)[c(19:20)] = c("external_gene_name_target_genes", "hgnc_symbol_target_genes")
merge_dbSUPER = merge_dbSUPER[c(1,3:12,2,13:22)]


# 4. TFs, diseases and mutations
merge_dbSUPER = cbind(merge_dbSUPER, external_gene_name_TFs = NA, hgnc_symbol_TFs = NA, TFs2enh_PMID = NA, TFs2enh_method = NA, 
                      disease = NA, disease_PMID = NA, disease_method = NA, refsnp_ID = NA, mutation_PMID = NA, mutation_method = NA)


#5. Filter
filt_merge_dbSUPER = merge_dbSUPER[-c(19,23)] #external gene name
filt_merge_dbSUPER = filt_merge_dbSUPER[!duplicated(filt_merge_dbSUPER),]

indexes = which(is.na(filt_merge_dbSUPER$current_chr))
if(length(indexes) > 0){filt_merge_dbSUPER = filt_merge_dbSUPER[-indexes,]} #sequences not mapped to hg38

indexes = which(is.na(filt_merge_dbSUPER$hgnc_symbol_target_gene)) #genes without gene symbol
if(length(indexes) > 0){filt_merge_dbSUPER$enh2gene_PMID[indexes] = filt_merge_dbSUPER$enh2gene_method[indexes] = NA}

filt_merge_dbSUPER[is.na(filt_merge_dbSUPER)] = "-"


#6. Write and save
rm(data, cell_name, index, indexes, PMID)
dir.create("./dbSUPER_results")
save.image("./dbSUPER_results/Data_dbSUPER.RData")
write.table(filt_merge_dbSUPER, "./dbSUPER_results/dbSUPER.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
