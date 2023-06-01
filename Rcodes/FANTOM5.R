#1. Reading files
library(data.table)
#Enhancers:
orig_enh_FANTOM5 = fread("https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz")
orig_enh_FANTOM5 = data.frame(orig_enh_FANTOM5)
enhancers = orig_enh_FANTOM5

#Enhancers-Target genes
orig_target_genes_FANTOM5 = fread("https://slidebase.binf.ku.dk/human_enhancers/presets/serve/enhancer_tss_associations", fill = T)
orig_target_genes_FANTOM5 = data.frame(orig_target_genes_FANTOM5)
target_genes_metadata = orig_target_genes_FANTOM5
target_genes_metadata[target_genes_metadata == ""] = NA
indexes = which(is.na(target_genes_metadata$V1))
if(length(indexes > 0)){target_genes_metadata = target_genes_metadata[-indexes,]}
#The file contained a blank line


#2. Enhancers
enhancers = enhancers[1:4]
colnames(enhancers) = c("orig_chr", "orig_start", "orig_end", "coord")
enhancer_metadata = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", enhancers$coord)), enhancers, orig_assembly = "hg19", 
                          original_ID = NA, crossref = NA, enh_PMID = 24670763, biosample_name = NA, enh_method = "CAGE", 
                          type = NA, source = "FANTOM5")
enh19 = enhancer_metadata[1:6]
rm(enhancers)

#hg38
source("./scripts/hgtohg_onestep.R")
hg19tohg38 = hgtohg_onestep(enh19[2:5], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
hg19tohg38_tables = hgtohg2tables("hg19", "hg38", enh19, hg19tohg38)

#Enhancer metadata 
enhancer_metadata = merge(hg19tohg38_tables$enhancers, enhancer_metadata[c(1,7:13)], by = "enh_ID", all = T)
enhancer_metadata = enhancer_metadata[!duplicated(enhancer_metadata),]


#3. Target genes
#Metadata
target_genes = target_genes_metadata$V4
coord = refseq = external_gene_name = R = FDR = c()
for (i in 1:length(target_genes)){
  split = unlist(strsplit(target_genes[i], ";"))
  coord[i] = split[1]
  refseq[i] = split[2]
  external_gene_name[i] = split[3]
  R[i] = split[4]
  FDR[i] = split[5]
}
target_genes_metadata = data.frame(TSS_chr = target_genes_metadata$V1, TSS_start = target_genes_metadata$V2, 
                                   TSS_end = target_genes_metadata$V3, coord, refseq, external_gene_name, R, FDR)
target_genes_metadata$R = unlist(lapply(target_genes_metadata$R, function(R){unlist(strsplit(R, ":"))[2]}))
target_genes_metadata$FDR = unlist(lapply(target_genes_metadata$FDR, function(FDR){unlist(strsplit(FDR, ":"))[2]}))

#Update gene symbol
target_genes = target_genes_metadata$external_gene_name
target_genes = unique(target_genes)
index = which(is.na(target_genes))
if(length(index) > 0){target_genes = target_genes[-index]}
source("./scripts/UpdateGeneSymbol.R")
target_genes = UpdateGeneSymbol(target_genes)
target_genes = target_genes[!duplicated(target_genes),]
target_genes_metadata = merge(target_genes_metadata, target_genes, by = "external_gene_name", all = T)
target_genes_metadata = target_genes_metadata[!duplicated(target_genes_metadata),]
target_genes_metadata = target_genes_metadata[c(2:5,1,9,6:8)]
source("./scripts/SplitIntoAtomicValues.R")
indexes = grep(",", target_genes_metadata$refseq)
data = target_genes_metadata[indexes,]
target_genes_metadata = target_genes_metadata[-indexes,]
data = SplitIntoAtomicValues(data, 7, ",")
target_genes_metadata = rbind(target_genes_metadata, data)
target_genes_metadata = target_genes_metadata[!duplicated(target_genes_metadata),]
target_genes_metadata = merge(target_genes_metadata, enh19[c(1,5)], by = "coord", all.x = T)
target_genes_metadata = target_genes_metadata[c(10,2:9)]
target_genes_metadata = cbind(target_genes_metadata, enh2gene_PMID = 24670763, 
                              enh2gene_method = "expression correlation") # within 500 kb window

FANTOM5 = merge(enhancer_metadata, target_genes_metadata[c(1,5:6,10:11)], by = "enh_ID", all = T)
FANTOM5 = FANTOM5[!duplicated(FANTOM5),]
colnames(FANTOM5)[19:20] = c("external_gene_name_target_genes", "hgnc_symbol_target_gene")


#4. TFs, diseases and mutations
FANTOM5 = cbind(FANTOM5, external_gene_name_TFs = NA, hgnc_symbol_TFs = NA, TFs2enh_PMID = NA, TFs2enh_method = NA, 
                disease = NA, disease_PMID = NA, disease_method = NA, refsnp_ID = NA, mutation_PMID = NA, mutation_method = NA)


#5. Filter
FANTOM5_filt = FANTOM5[-c(19, 23)]
FANTOM5_filt = FANTOM5_filt[!duplicated(FANTOM5_filt),]

indexes = which(is.na(FANTOM5_filt$current_chr))
if(length(indexes) > 0){FANTOM5_filt = FANTOM5_filt[-indexes,]}

indexes = which(is.na(FANTOM5_filt$hgnc_symbol_target_gene))
if(length(indexes) > 0){FANTOM5_filt$enh2gene_PMID[indexes] = FANTOM5_filt$enh2gene_method[indexes] = NA}

FANTOM5_filt[is.na(FANTOM5_filt)] = "-"


#6. Write and save
rm(data, coord, external_gene_name, FDR, i, index, indexes, R, refseq, split)
dir.create("./FANTOM5_results")
save.image("./FANTOM5_results/Data_FANTOM5.RData")
write.table(FANTOM5_filt, "./FANTOM5_results/FANTOM5.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
