#1. Read files
#ENCODE+Roadmap Lasso
#http://yiplab.cse.cuhk.edu.hk/jeme/encoderoadmap_lasso.zip
library(data.table)
encoderoadmap_lasso = data.frame()
path = "./JEME_data/encoderoadmap_lasso/"
list_files = list.files(path)
for (i in 1:length(list_files)){
  biosample = unlist(strsplit(list_files[i], "\\."))[2]
  data = fread(paste0(path, list_files[i]))
  data = cbind(data, biosample)
  encoderoadmap_lasso = rbind(encoderoadmap_lasso, data)
}
rm(data, biosample, i, list_files, path)

#ENCODE+Roadmap Elastic Net
#http://yiplab.cse.cuhk.edu.hk/jeme/encoderoadmap_elasticnet.zip
encoderoadmap_elasticnet = data.frame()
path = "./JEME_data/encoderoadmap_elasticnet/"
list_files = list.files(path)
for (i in 1:length(list_files)){
  biosample = unlist(strsplit(list_files[i], "\\."))[2]
  data = fread(paste0(path, list_files[i]))
  data = cbind(data, biosample)
  encoderoadmap_elasticnet = rbind(encoderoadmap_elasticnet, data)
}
rm(data, biosample, i, list_files, path)

#FANTOM5 Lasso
#http://yiplab.cse.cuhk.edu.hk/jeme/fantom5_lasso.zip
fantom5_lasso = data.frame()
path = "./JEME_data/fantom5_lasso/"
list_files = list.files(path)
for (i in 1:length(list_files)){
  biosample = unlist(strsplit(list_files[i], "\\."))[2]
  data = fread(paste0(path, list_files[i]))
  data = cbind(data, biosample)
  fantom5_lasso = rbind(fantom5_lasso, data)
}
rm(data, biosample, i, list_files, path)

#FANTOM5 Elastic Net
#http://yiplab.cse.cuhk.edu.hk/jeme/fantom5_elasticnet.zip
fantom5_elasticnet = data.frame()
path = "./JEME_data/fantom5_elasticnet/"
list_files = list.files(path)
for (i in 1:length(list_files)){
  biosample = unlist(strsplit(list_files[i], "\\."))[2]
  data = fread(paste0(path, list_files[i]))
  data = cbind(data, biosample)
  fantom5_elasticnet = rbind(fantom5_elasticnet, data)
}
rm(data, biosample, i, list_files, path)


#2. Enh19 to enh38
enh19 = c(encoderoadmap_lasso$V1, encoderoadmap_elasticnet$V1, fantom5_lasso$V1, fantom5_elasticnet$V1)
enh19 = unique(enh19)
source("./scripts/coord2bed.R")
enh19 = coord2bed(enh19)
enh_ID = paste0("hg19_", gsub("[:-]", "_", enh19$coord))
enh19 = cbind(enh_ID, enh19, orig_assembly = "hg19")
colnames(enh19)[2:4] = c("orig_chr", "orig_start", "orig_end")

source("./scripts/hgtohg_onestep.R")
hg19tohg38 = hgtohg_onestep(enh19[2:6], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
hg19tohg38_tables = hgtohg2tables("hg19", "hg38", enh19, hg19tohg38)


#3.Enhancer datatable
encoderoadmap_lasso = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", encoderoadmap_lasso$V1)), encoderoadmap_lasso, 
                                            enh2gene_method = "Lasso", enh_method = "ChIP-seq,ChromHMM")
encoderoadmap_elasticnet = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", encoderoadmap_elasticnet$V1)), 
                                 encoderoadmap_elasticnet, enh2gene_method = "Elastic Net", enh_method = "ChIP-seq,ChromHMM")
tables = rbind(encoderoadmap_lasso, encoderoadmap_elasticnet)

library(readxl)
encoderoadmap_samples = read_excel("./JEME_data/biosamples.xlsx", sheet = "encoderoadmap")
colnames(encoderoadmap_samples)[1] = "biosample"
encoderoadmap_samples$biosample = as.character(encoderoadmap_samples$biosample)
enh_JEME = merge(tables, encoderoadmap_samples[,c(1,5)], by = "biosample", all.x = T)
enh_JEME = enh_JEME[,-1]
colnames(enh_JEME)[7] = "biosample_name"

fantom5_lasso = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", fantom5_lasso$V1)), fantom5_lasso, enh2gene_method = "Lasso",
                      enh_method = "CAGE")
fantom5_elasticnet = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", fantom5_elasticnet$V1)), fantom5_elasticnet, 
                           enh2gene_method = "Elastic Net", enh_method = "CAGE")
tables = rbind(fantom5_lasso, fantom5_elasticnet)
fantom5_samples = read_excel("./JEME_data/biosamples.xlsx", sheet = "fantom5")
colnames(fantom5_samples)[1] = "biosample"
fantom5_samples$biosample = as.character(fantom5_samples$biosample)
tables = merge(tables, fantom5_samples[,c(1,3)], by = "biosample", all.x = T)
tables = tables[,-1]
colnames(tables)[7] = "biosample_name"

enh_JEME = rbind(enh_JEME, tables)
enh_JEME = merge(hg19tohg38_tables$enhancers, enh_JEME, by = "enh_ID", all.x = T)
rm(tables, enh_ID)
enh_JEME = enh_JEME[!duplicated(enh_JEME),]

colnames(enh_JEME)[12:14] = c("coord", "external_target_gene_names", "score_tgene")
enh_JEME = data.frame(enh_JEME[,1:11], original_ID = NA, crossref = NA, enh_PMID = 28869592, biosample_name = enh_JEME$biosample_name,
                      enh_method = enh_JEME$enh_method, type = NA, source = "JEME", external_target_gene_names = enh_JEME$external_target_gene_names,
                      enh2gene_PMID = 28869592, enh2gene_method = enh_JEME$enh2gene_method, hgnc_symbol_TFs = NA, TFs2enh_PMID = NA,
                      TFs2enh_method = NA, disease = NA, disease_PMID = NA, disease_method = NA, refseq_ID = NA, mutation_PMID = NA, 
                      mutation_method = NA)
indexes = grep(",", enh_JEME$enh_method)
set = enh_JEME[indexes,]
enh_JEME$enh_method[indexes] = "ChIP-seq"
set$enh_method = "ChromHMM"
enh_JEME = rbind(enh_JEME, set)
rm(set, indexes)


#4. Updating genes
external_target_gene_names = unlist(lapply(enh_JEME$external_target_gene_names, function(gene){unlist(strsplit(gene, "\\$"))[2]}))
table_genes = as.data.frame(cbind(enh_JEME$external_target_gene_names, external_target_gene_names))
rm(external_target_gene_names)
table_genes = table_genes[!duplicated(table_genes),]
colnames(table_genes) = c("external_target_gene_names", "external_gene_name")

external_gene_names = unique(table_genes$external_gene_name)
source("./scripts/UpdateGenes.R")
upd_genes = UpdateGenes(external_gene_names)
table_genes = merge(table_genes, upd_genes, by = "external_gene_name", all.x = T)
table_genes = table_genes[!duplicated(table_genes),]

enh_JEME = merge(enh_JEME, table_genes[,2:3], by = "external_target_gene_names", all.x = T)
colnames(enh_JEME)[31] = "hgnc_symbol_target_genes"
enh_JEME = enh_JEME[,c(2:19,1,31,20:30)]
enh_JEME = enh_JEME[!duplicated(enh_JEME),]


#5. Filtering data
filt_JEME = enh_JEME[,-19]
indexes = which(is.na(filt_JEME$hgnc_symbol_target_genes))
filt_JEME$enh2gene_PMID[indexes] = filt_JEME$enh2gene_method[indexes] = NA
indexes = which(is.na(filt_JEME$current_chr))
filt_JEME = filt_JEME[-indexes,]
filt_JEME[is.na(filt_JEME)] = "-"
rm(indexes)


#5. Saving files
dir.create("./JEME_results/")
save.image("./JEME_results/JEME_data.RData")
write.table(filt_JEME, "./JEME_results/JEME.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
write(sort(unique(filt_JEME$biosample_name)), "./biosamples.tsv", sep = "\n")

