#1. Load files: http://www.enhanceratlas.org/downloadv2.ph
library(data.table)
list_files = list.files("./EnhancerAtlas2.0_data/Enhancer/")
EnhancerAtlas = data.frame()
for (i in 1:length(list_files)){
  biosample = unlist(strsplit(list_files[i], "\\."))[1]
  path = paste0("./EnhancerAtlas2.0_data/Enhancer/", list_files[i])
  file = fread(path)
  original_ID = paste0("HS", paste(rep(0, 3-nchar(i)), collapse = ""), i, "-", sapply(seq(1:nrow(file)), function(x){
    paste0(paste(rep(0,5-nchar(x)), collapse = ""), x)
  }))
  file = cbind(file, biosample_name = biosample, original_ID)
  EnhancerAtlas = rbind(EnhancerAtlas, file)
}
rm(file, biosample, i, list_files, path, original_ID)
EnhancerAtlas = EnhancerAtlas[!duplicated(EnhancerAtlas),]


#2. Enhancer Table
source("./scripts/bed2coord.R")
coord = bed2coord(EnhancerAtlas[,1:3])
EnhancerAtlas = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", coord)), EnhancerAtlas[,1:3], 
                      coord, orig_assembly = "hg19", original_ID = EnhancerAtlas$original_ID, 
                      crossref = paste0("http://www.enhanceratlas.org/browseenhancer.php?enhancer=", EnhancerAtlas$original_ID), 
                      enh_PMID = 31740966, biosample_name = EnhancerAtlas$biosample_name, 
                      enh_method = "unsupervised machine learning", type = "-", 
                      source = "EnhancerAtlas2.0", orig_score = EnhancerAtlas$V4)
colnames(EnhancerAtlas)[2:4] = c("orig_chr", "orig_start", "orig_end")
rm(coord)

#hg19 to hg38
source("./scripts/hgtohgparallel.R")
library(parallel)
enh19 = EnhancerAtlas[,1:6]
enh19 = enh19[!duplicated(enh19),]
hg19tohg38 = hgtohgparallel(enh19[2:5], 19, 38, detectCores()/2)

source("./scripts/size_comparison_parallel.R")
size_hg2hg = size_comparison_parallel(hg19tohg38, 0.05)

source("./scripts/filtAfterSizeComp_parallel.R")
filtered = filtAfterSizeComp_parallel(hg19tohg38, size_hg2hg, enh19[2:5])
enhancers = cbind(filtered$merge_table, score = size_hg2hg$score)
indexes = which(is.na(enhancers$new_coord))
enhancers[indexes, "score"] = NA
rm(indexes)

enh_ID = paste0("hg19_", gsub("[:-]", "_", enhancers$orig_coord))
enhancers = cbind(enh_ID, enhancers)
EnhancerAtlas = merge(EnhancerAtlas, enhancers[,c(1,6:8,10)], all.x = T)
EnhancerAtlas = cbind(EnhancerAtlas[,c(1:4,6,15:17)], current_assembly = "hg38", minimum_ratio = 0.95, EnhancerAtlas[,c(18,7:14)])
colnames(EnhancerAtlas)[6:8] = c("current_chr", "current_start", "current_end")


#3. Interactions
list_files = list.files("./EnhancerAtlas2.0_data/Interactions/")
EnhancerAtlas_EP = data.frame()
for (i in 1:length(list_files)){
  biosample = unlist(strsplit(list_files[i], "\\."))[1]
  biosample = unlist(strsplit(biosample, "_"))[1]
  path = paste0("./EnhancerAtlas2.0_data/Interactions/", list_files[i])
  file = fread(path)
  split = unlist(strsplit(file$V1, "\\$"))
  file = cbind(data.frame(matrix(split, nrow = nrow(file), byrow = T), stringsAsFactors = F), file[,2])
  split = unlist(strsplit(file$X1, "_"))
  file = cbind(data.frame(matrix(split, nrow = nrow(file), byrow = T), stringsAsFactors = F), file[,2:ncol(file)],
               biosample_name = biosample)
  EnhancerAtlas_EP = rbind(EnhancerAtlas_EP, file)
}
colnames(EnhancerAtlas_EP) = c("enh", "ensembl_gene", "external_target_gene_name", "chr", "TSS", "strand", "score", "biosample_name")
EnhancerAtlas_EP = cbind(EnhancerAtlas_EP, enh2gene_PMID = 31740966, enh2gene_method = "EAGLE")
rm(list_files, i, biosample, path, file, split)


#Update genes
genes = unique(EnhancerAtlas_EP$gene_symbol)
source("./scripts/UpdateGene.R")
upd_genes = UpdateGeneSymbol(genes)
colnames(upd_genes) = c("external_target_gene_name", "hgnc_symbol_target_genes")
EnhancerAtlas_EP = merge(EnhancerAtlas_EP, upd_genes, by = "external_target_gene_name", all.x = T)
EnhancerAtlas_EP = EnhancerAtlas_EP[!duplicated(EnhancerAtlas_EP),]
EnhancerAtlas_EP = EnhancerAtlas_EP[,c(2:8,1,11,9,10)]
EnhancerAtlas_EP = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", EnhancerAtlas_EP$enh)), EnhancerAtlas_EP)


#4. Complete set
enh_EnhancerAtlas = merge(EnhancerAtlas[,-19], EnhancerAtlas_EP[,c(1,8:12)], by = c("enh_ID", "biosample_name"), all.x = T)
enh_EnhancerAtlas = enh_EnhancerAtlas[!duplicated(enh_EnhancerAtlas),]
enh_EnhancerAtlas = cbind(enh_EnhancerAtlas[,c(1,3:15,2,16:22)], hgnc_symbol_TFs = NA, TFs2enh_PMID = NA, TFs2enh_method = NA, 
                          disease = NA, disease_PMID = NA, disease_method = NA, 
                          refseq_ID = NA, mutation_PMID = NA, mutation_method = NA)


#5. Filtered table
filt_EnhancerAtlas = enh_EnhancerAtlas[,-19]
indexes = which(is.na(filt_EnhancerAtlas$hgnc_symbol_target_genes))
filt_EnhancerAtlas$enh2gene_method[indexes] = filt_EnhancerAtlas$enh2gene_PMID[indexes] = NA
indexes = which(is.na(filt_EnhancerAtlas$current_chr))
filt_EnhancerAtlas = filt_EnhancerAtlas[-indexes,]
filt_EnhancerAtlas$enh2gene_method = as.character(filt_EnhancerAtlas$enh2gene_method)
filt_EnhancerAtlas[is.na(filt_EnhancerAtlas)] = "-"
rm(indexes)


#6. Save files
dir.create("./EnhancerAtlas2.0_results/")
save.image("./EnhancerAtlas2.0_results/EnhancerAtlas2.0_data.RData")
write.table(filt_EnhancerAtlas, "./EnhancerAtlas2.0_results/EnhancerAtlas.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
