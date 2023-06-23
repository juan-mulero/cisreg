#1. Enhancers

library(data.table)
options(scipen=99999)

# http://enhanceratlas.net/scenhancer/data/download/enhancer/hs/

path = "./scEnhancer_data/enhancers/"
enh_files = list.files(path)
enhancers = data.frame()
for (i in 1:length(enh_files)){
  data = fread(paste0(path, enh_files[i]), header = F)
  colnames(data) = c("orig_chr", "orig_start", "orig_end", "orig_score")
  biosample = gsub(".bed", "", enh_files[i])
  data = cbind(data, biosample_name = biosample)
  enhancers = rbind(enhancers, data)
}
rm(data, biosample, enh_files, i, path)

source("./scripts/bed2coord.R")
coord = bed2coord(enhancers[,1:3])
enh_ID = paste0("hg19_", gsub("[:-]", "_", coord))
enhancers = cbind(enh_ID, enhancers[,1:3], orig_assembly = "hg19", orig_score = enhancers$orig_score, original_ID = "-",
                  crossref = "-", enh_PMID = 34761274, biosample_name = enhancers$biosample_name, 
                  enh_method = "ATAC-seq,unsupervised machine learning", type = "-", source = "scEnhancer")
enhancers = enhancers[!duplicated(enhancers),]
rm(enh_ID, coord)


#2. Enh-Target gene

# http://enhanceratlas.net/scenhancer/data/download/TargetGene/human/

path = "./scEnhancer_data/TargetGene/"
target_genes_files = list.files(path)
target_genes = data.frame()
for (i in 1:length(target_genes_files)){
  data = fread(paste0(path, target_genes_files[i]), header = F)
  colnames(data) = c("enh_coord", "target_gene")
  biosample = gsub(".bed", "", target_genes_files[i])
  data = cbind(data, biosample_name = biosample)
  target_genes = rbind(target_genes, data)
}
rm(data, biosample, target_genes_files, i, path)

enh_ID = paste0("hg19_", gsub("[:-]", "_", target_genes$enh_coord))
source("./scripts/coord2bed.R")
bed_table = coord2bed(target_genes$enh_coord)
colnames(bed_table) = c("orig_chr", "orig_start", "orig_end", "coord")
target_genes = cbind(enh_ID = enh_ID, bed_table[,-4], orig_assembly = "hg19", target_genes[,-1], 
                     enh2gene_PMID = 34761274, enh2gene_method = "Cicero")
rm(enh_ID, bed_table)

source("./scripts/SplitColumn.R")
target_genes = SplitColumn(target_genes, 6, ";")

#We split the column target_gene into 2 column
new_target_gene = strsplit(target_genes$target_gene, ":")
df_new_target_genes = data.frame()
n = 50000
n_groups = ceiling(length(new_target_gene)/n)
init = 1
t = proc.time()
for (group in 1:n_groups){
  cat("Group", group, "of", n_groups, "\n", proc.time() - t, "\n")
  if (group == n_groups){
    final = length(new_target_gene)
  } else {
    final = init + n - 1
  }
  set = t(data.frame(new_target_gene[init:final]))
  df_new_target_genes = rbind(df_new_target_genes, set)
  init = final + 1
}
colnames(df_new_target_genes) = c("ensembl_gene_id", "external_gene_name")
target_genes = cbind(target_genes, df_new_target_genes)
rm(df_new_target_genes, new_target_gene, set, final, group, init, n, n_groups, t)
target_genes = target_genes[,-6]


#3. hg19 to hg38 and complete set 
scEnhancers = merge(enhancers, target_genes, by = c("enh_ID", "orig_chr", "orig_start", "orig_end", 
                                                    "orig_assembly", "biosample_name"), all = T)
scEnhancers$enh_PMID = 34761274
scEnhancers$enh_method = "ATAC-seq,unsupervised machine learning"
scEnhancers$source = "scEnhancer"
scEnhancers[is.na(scEnhancers)] = "-"
scEnhancers = scEnhancers[!duplicated(scEnhancers),]

coord = bed2coord(scEnhancers[,2:4])
bed_table = data.frame(scEnhancers[,2:4], coord)
bed_table = bed_table[!duplicated(bed_table),]
colnames(bed_table)[4] = "coord"
bed_table = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", bed_table$coord)), bed_table, assembly = "hg19")
colnames(bed_table) = c("enh_ID", "orig_chr", "orig_start", "orig_end", "coord", "orig_assembly")
rm(coord, enhancers)

source("./scripts/hgtohg_onestep.R")
hg19tohg38 = hgtohg_onestep(bed_table[,2:5], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
hg19tohg38_tables = hgtohg2tables("hg19", "hg38", bed_table, hg19tohg38)
scEnhancers = merge(hg19tohg38_tables$enhancers, scEnhancers[,c(1,6:17)], by = "enh_ID", all = T)
scEnhancers = scEnhancers[!duplicated(scEnhancers),]
scEnhancers[is.na(scEnhancers)] = "-"
rm(bed_table)


#4. SNPs

# http://enhanceratlas.net/scenhancer/data/download/SNP/hs/

path = "./scEnhancer_data/SNPs/"
snp_files = list.files(path)
snps = data.frame()
for (i in 1:length(snp_files)){
  data = fread(paste0(path, snp_files[i]), header = F)
  colnames(data) = c("orig_chr", "orig_start", "orig_end", "orig_score", "snp_chr", "snp_start", "snp_end", "refsnp_ID")
  biosample = gsub(".bed", "", snp_files[i])
  data = cbind(data, biosample_name = biosample)
  snps = rbind(snps, data)
}
rm(data, biosample, snp_files, i, path)

coord = bed2coord(snps[,1:3])
snps = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", coord)), snps)
rm(coord)
scEnhancers = merge(scEnhancers, snps[,c(1,5:10)], by = c("enh_ID", "biosample_name", "orig_score"), all = T)
scEnhancers[is.na(scEnhancers)] = "-"
scEnhancers = scEnhancers[!duplicated(scEnhancers),]
scEnhancers = scEnhancers[,c(1,4:16,2,17,3,18,19,22,23,20,21,27,24:26)]
scEnhancers = cbind(scEnhancers[,1:23], hgnc_symbol_TFs = "-", TFs2enh_PMID = "-", TFs2enh_method = "-",
                    disease = "-", disease_PMID = "-", disease_method = "-", scEnhancers[,24:27])


#5. Interactions

# http://enhanceratlas.net/scenhancer/data/download/interaction/hs/

path = "./scEnhancer_data/Interactions/"
int_files = list.files(path)
interactions = data.frame()
for (i in 1:length(int_files)){
  data = fread(paste0(path, int_files[i]), header = F, sep = "\t")
  colnames(data) = c("orig_enh", "orig_target", "int_score")
  biosample = gsub("_interaction.txt", "", int_files[i])
  data = cbind(data, biosample_name = biosample)
  interactions = rbind(interactions, data)
}
rm(data, biosample, int_files, i, path)

orig_enh = strsplit(interactions$orig_enh, "\\|")
orig_enh = t(data.frame(orig_enh))

orig_target = strsplit(interactions$orig_target, "\\|")
orig_target = t(data.frame(orig_target))

interactions = cbind(orig_enh, orig_target, interactions[,3:4])
rm(orig_enh, orig_target)
colnames(interactions) = c("orig_enh", "orig_score", "orig_target", "target_gene", "enh_gene_score", "biosample_name")
interactions$orig_enh = paste0("hg19_", gsub("[:-]", "_", interactions$orig_enh))
colnames(interactions)[1] = "enh_ID"

target_gene = unlist(lapply(interactions$target_gene, function(gene){unlist(strsplit(gene, ":"))[1]}))
interactions$target_gene = target_gene
colnames(interactions)[4] = "external_target_gene_name"
colnames(scEnhancers)[21] = "external_target_gene_name"
rm(target_gene)

scEnhancers = merge(scEnhancers, interactions[,-2], by = c("enh_ID", "external_target_gene_name", "biosample_name"),
               all.x = T)
scEnhancers = scEnhancers[!duplicated(scEnhancers),]
scEnhancers[is.na(scEnhancers)] = "-"
scEnhancers = scEnhancers[,c(1,4:16,3,17:20,34,2,21:23,35,24:33)]
scEnhancers = cbind(scEnhancers, mutation_PMID = "-", mutation_method = "-")
indexes = which(scEnhancers$refsnp_ID != "-")
scEnhancers$mutation_PMID[indexes] = 34761274
rm(indexes)


#6. Updating genes
genes = unique(scEnhancers$external_target_gene_name)
index = which(genes == "-")
genes = genes[-index]
source("./scripts/UpdateGenes.R")
upd_genes = UpdateGenes(genes)
colnames(upd_genes) = c("external_target_gene_name", "hgnc_symbol_target_genes")
rm(genes, index)
scEnhancers = merge(scEnhancers, upd_genes, by = "external_target_gene_name", all.x = T)
scEnhancers[is.na(scEnhancers)] = "-"
scEnhancers = scEnhancers[!duplicated(scEnhancers),]
scEnhancers = scEnhancers[,c(2:21,1,38,22:37)]
scEnhancers = scEnhancers[order(scEnhancers$enh_ID),]


# 7. Filtering and saving
indexes = which(scEnhancers$current_chr == "-")
filt_scEnhancers = scEnhancers[-indexes,-c(17,20,21,23,26,34:36)] 
indexes = which(filt_scEnhancers$hgnc_symbol_target_genes == "-")
filt_scEnhancers$enh2gene_PMID[indexes] = filt_scEnhancers$enh2gene_method[indexes] = "-"
rm(indexes)

set = filt_scEnhancers
set$enh_method = "unsupervised machine learning"
filt_scEnhancers$enh_method = "ATAC-seq"
filt_scEnhancers = rbind(filt_scEnhancers, set)
rm(set)
filt_scEnhancers = filt_scEnhancers[order(filt_scEnhancers$enh_ID),]

dir.create("./scEnhancer_results/")
save.image("./scEnhancer_results/scEnhancer_data.RData")
write.table(filt_scEnhancers, "./scEnhancer_results/scEnhancer.tsv", col.names = F, row.names = F, quote = F, sep = "\t")
