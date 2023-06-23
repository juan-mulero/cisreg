#1. Reading file
library(data.table)
orig_HACER = data.frame(fread("http://bioinfo.vanderbilt.edu/AE/HACER/download/T1.txt", header = F))
colnames(orig_HACER) = c("original_ID", "orig_chr", "orig_start", "orig_end", "orig_center", "IsInFANTOM5", "geneFANTOM5", "gene50kb",
                         "gene4DGenome", "biosample_name_4DGenome", "method_4DGenome", "PMID_4DGenome", "closest_gene", 
                         "distance_closest_gene", "enh_method", "biosample_name", "orig_assembly", "source", "norm_enh_expr", 
                         "log2_density", "IsInVISTA", "IsInEnsembl", "IsInENCODE", "IsInChromHMM")
source("./scripts/bed2coord.R")
coord = bed2coord(orig_HACER[2:4])
enh_ID = paste0("hg19_", gsub("[:-]", "_", coord))
crossref = paste0("http://bioinfo.vanderbilt.edu/AE/HACER/detail.php?id=", orig_HACER$original_ID)
enh_metadata = cbind(enh_ID, original_ID = orig_HACER$original_ID, crossref, enh_PMID = 30247654, 
                     biosample_name = orig_HACER$biosample_name, enh_method = orig_HACER$enh_method, type = NA, source = "HACER")
enh_metadata = data.frame(enh_metadata)


#2. hg19 to hg38
source("./scripts/hgtohg_onestep.R")
enh19 = cbind(enh_ID, orig_HACER[c(2:4)], coord, orig_HACER[17])
enh19 = enh19[!duplicated(enh19),]
enh19to38 = hgtohg_onestep(enh19[2:5], 19, 38, 0.05)
source("./scripts/hgtohg2tables.R")
tables_HACER = hgtohg2tables("hg19", "hg38", enh19, enh19to38)
enh_HACER = tables_HACER$enhancers
enh_HACER = merge(enh_HACER, enh_metadata, by = "enh_ID")
rm(coord, crossref, enh_ID, enh19, enh_metadata)


#3. Genes
source("./scripts/SplitColumn.R")
split_orig_HACER = orig_HACER
split_orig_HACER$geneFANTOM5 = gsub("NA;|;NA", "", split_orig_HACER$geneFANTOM5)
split_orig_HACER$geneFANTOM5 = gsub(";", ",", split_orig_HACER$geneFANTOM5)

split_orig_HACER = SplitColumn(split_orig_HACER, 7, ",") #geneFANTOM5
split_orig_HACER = split_orig_HACER[!duplicated(split_orig_HACER),]
split_orig_HACER = SplitColumn(split_orig_HACER, 8, ",") #gene50kb
split_orig_HACER = split_orig_HACER[!duplicated(split_orig_HACER),]
split_orig_HACER = SplitColumn(split_orig_HACER, 13, ",") #closest
split_orig_HACER = split_orig_HACER[!duplicated(split_orig_HACER),]

external_gene_name_target_genes = c(split_orig_HACER$geneFANTOM5, split_orig_HACER$gene50kb, split_orig_HACER$closest_gene)
external_gene_name_target_genes = unique(external_gene_name_target_genes)
indexes = which(is.na(external_gene_name_target_genes) | external_gene_name_target_genes == "NA")
if (length(indexes) > 0){external_gene_name_target_genes = external_gene_name_target_genes[-indexes]}
source("./scripts/UpdateGeneSymbol.R")
upd_target_genes = UpdateGeneSymbol(external_gene_name_target_genes)
colnames(upd_target_genes) = c("external_gene_name_target_genes", "hgnc_symbol_target_genes")

#geneFANTOM5
set1 = cbind(split_orig_HACER[c(1,7)], enh2gene_method = "expression correlation", enh2gene_PMID = 30247654)
#gene50kb
set2 = cbind(split_orig_HACER[c(1,8)], enh2gene_method = "distance_50kb", enh2gene_PMID = 30247654)
#closest_gene
set3 = cbind(split_orig_HACER[c(1,13)], enh2gene_method = "closest", enh2gene_PMID = 30247654)

colnames(set1)[2] = colnames(set2)[2] = colnames(set3)[2] = "external_gene_name_target_genes"
set = rbind(set1, set2, set3)
rm(set1, set2, set3)
indexes = which(is.na(set$external_gene_name_target_genes))
set = set[-indexes,]
set = merge(set, upd_target_genes, by = "external_gene_name_target_genes", all.x = T)
set = set[!duplicated(set),]
enh_HACER = merge(enh_HACER, set, by = "original_ID", all.x = T)
enh_HACER = enh_HACER[!duplicated(enh_HACER),]
enh_HACER = cbind(enh_HACER[c(2:12,1,13:19,22,21,20)], hgnc_symbol_TFs = NA, TFs2enh_PMID = NA, TFs2method = NA, disease = NA,
                  disease_PMID = NA, disease_method = NA, refsnp_ID = NA, mutation_PMID = NA, mutation_method = NA)
rm(set, indexes)


#4. Filtered set
filt_HACER = enh_HACER[,-19]
indexes = which(is.na(filt_HACER$hgnc_symbol_target_genes))
filt_HACER$enh2gene_method[indexes] = filt_HACER$enh2gene_PMID[indexes] = NA
indexes = which(is.na(filt_HACER$current_chr))
filt_HACER = filt_HACER[-indexes,]
filt_HACER = data.frame(filt_HACER)
filt_HACER$enh2gene_method = as.character(filt_HACER$enh2gene_method)
filt_HACER[is.na(filt_HACER)] = "-"

indexes = which(filt_HACER$orig_chr == "chrM")
filt_HACER = filt_HACER[-indexes,]
rm(indexes)


#5. Save file
dir.create("./HACER_results")
save.image("./HACER_results/data_HACER.RData")
write.table(filt_HACER, "./HACER_results/HACER.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
