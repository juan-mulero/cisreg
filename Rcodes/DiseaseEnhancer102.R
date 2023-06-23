#1. Data collection
library(data.table)
disease_enh = fread("http://biocc.hrbmu.edu.cn/DiseaseEnhancer/RFunctions/enh2disease-1.0.2.txt")
colnames(disease_enh) = c("id", "chr", "start", "end", "TargetGene", "DiseaseType")

mutations = fread("http://biocc.hrbmu.edu.cn/DiseaseEnhancer/RFunctions/enhInfo-1.0.2.txt")
colnames(mutations) = c("id", "VariantType", "VariantName", "chr", "start", "end", "TargetGene", "MutationType", "VariantConsequence", "PMID")

enh_target = fread("http://biocc.hrbmu.edu.cn/DiseaseEnhancer/RFunctions/enh2target-1.0.2.txt")
colnames(enh_target) = c("chr", "start", "end", "TargetGene")


#2. Enhancer table
disease_enh$chr[disease_enh$chr == "chrx"] = "chrX"

#hg19
source("./scripts/bed2coord.R")
coord = bed2coord(disease_enh[,c(2:4)])
enh_ID = paste0("hg19_", gsub("[:-]", "_", coord))
hg19 = cbind(enh_ID, disease_enh[,c(2:4)], coord, orig_assembly = "hg19")
hg19 = hg19[!duplicated(hg19),]
colnames(hg19)[2:4] = c("orig_chr", "orig_start", "orig_end")

#hg19 to hg38
source("./scripts/hgtohg_onestep.R")
indexes = which(is.na(hg19$orig_chr))
hg19 = hg19[-indexes,]
hg19tohg38 = hgtohg_onestep(hg19[,2:5], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
tables_DiseaseEnhancer = hgtohg2tables("hg19", "hg38", hg19, hg19tohg38)
enh_DiseaseEnhancer = tables_DiseaseEnhancer$enhancers


#3. Enrich table
#Enrich genes
disease_enh = cbind(enh_ID, disease_enh)
enh_DiseaseEnhancer = merge(enh_DiseaseEnhancer, disease_enh[,c(1,2,6,7)], by = "enh_ID", all.x = T)
enh_DiseaseEnhancer = enh_DiseaseEnhancer[!duplicated(enh_DiseaseEnhancer),]
colnames(enh_DiseaseEnhancer)[12:14] = c("original_ID", "external_target_gene_name", "disease_name")

source("./scripts/SplitColumn.R")
enh_DiseaseEnhancer = SplitColumn(enh_DiseaseEnhancer, 13, ", ")

#Enrich PMID and mutations
mutations = merge(mutations, disease_enh[,1:2], by = "id", all.x = T)
mutations = mutations[!duplicated(mutations),]
enh_DiseaseEnhancer = merge(data.frame(enh_DiseaseEnhancer), data.frame(mutations[,c(3,10:11)]), by = "enh_ID", all.x = T)
enh_DiseaseEnhancer = enh_DiseaseEnhancer[!duplicated(enh_DiseaseEnhancer),]
enh_DiseaseEnhancer = SplitColumn(enh_DiseaseEnhancer, 16, ", ")
enh_DiseaseEnhancer = enh_DiseaseEnhancer[!duplicated(enh_DiseaseEnhancer),]

enh_DiseaseEnhancer = cbind(enh_DiseaseEnhancer[,1:12], crossref = NA, enh_PMID = enh_DiseaseEnhancer$PMID, biosample_name = NA,
                            enh_method = NA, type = NA, source = "DiseaseEnhancer", external_target_gene_name = enh_DiseaseEnhancer[,13],
                            enh2gene_PMID = enh_DiseaseEnhancer$PMID, enh2gene_method = NA, hgnc_symbol_TFs = NA,
                            TFs2enh_PMID = NA, TFs2enh_method = NA, disease_name = enh_DiseaseEnhancer[,14], 
                            disease_PMID = enh_DiseaseEnhancer$PMID, disease_method = NA, refseq_ID = enh_DiseaseEnhancer$VariantName, 
                            mutation_PMID = enh_DiseaseEnhancer$PMID, mutation_method = NA)
enh_DiseaseEnhancer[enh_DiseaseEnhancer == ""] = NA
indexes = which(grepl("^rs", enh_DiseaseEnhancer$refseq_ID) == F)
enh_DiseaseEnhancer$refseq_ID[indexes] = enh_DiseaseEnhancer$mutation_PMID[indexes] = NA
rm(coord, enh_ID, indexes)


#4. Genes
source("./scripts/UpdateGeneSymbol.R")
enh_DiseaseEnhancer$external_target_gene_name = gsub(" ", "-", enh_DiseaseEnhancer$external_target_gene_name)
genes = unique(enh_DiseaseEnhancer$external_target_gene_name)
upd_genes = UpdateGeneSymbol(genes)
colnames(upd_genes) = c("external_target_gene_name", "hgnc_symbol_target_genes")
enh_DiseaseEnhancer = merge(enh_DiseaseEnhancer, upd_genes, by = "external_target_gene_name", all.x = T)


#5. Filter data
enh_DiseaseEnhancer = enh_DiseaseEnhancer[,c(2:19,1,31,20:30)]
indexes = which(is.na(enh_DiseaseEnhancer$hgnc_symbol_target_genes))
filt_DiseaseEnhancer = enh_DiseaseEnhancer[,-19]
filt_DiseaseEnhancer$enh2gene_PMID[indexes] = NA
filt_DiseaseEnhancer[is.na(filt_DiseaseEnhancer)] = "-"
rm(indexes)


#6. Disease mapping
library(readxl)
library(tidyr)
diseases = read_excel("diseases.xlsx")
diseases = pivot_longer(diseases[1:4], cols = c("OMIM", "DiseaseOntology", "MeSH"), names_to = "ontology", values_to = "ID")
diseases = diseases[!duplicated(diseases),]
indexes = which(diseases$ID == "-")
diseases = diseases[-indexes, -2]
colnames(diseases)[1] = "disease_name"
filt_DiseaseEnhancer = merge(filt_DiseaseEnhancer, diseases, by = "disease_name", all.x = T)
filt_DiseaseEnhancer = filt_DiseaseEnhancer[,c(2:25,31,26:30)]
colnames(filt_DiseaseEnhancer)[25] = "disease"
filt_DiseaseEnhancer = filt_DiseaseEnhancer[!duplicated(filt_DiseaseEnhancer),]
indexes = which(filt_DiseaseEnhancer$disease == "-")
rm(indexes)


#6. Save files
dir.create("./DiseaseEnhancer102_results/")
save.image("./DiseaseEnhancer102_results/DiseaseEnhancer102_data.RData")
write.table(filt_DiseaseEnhancer, "./DiseaseEnhancer102_results/DiseaseEnhancer.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
