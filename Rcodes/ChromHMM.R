#1. Read files
biosamples = c("GM12878", "H1-hESC", "HepG2", "HMEC", "HSMM", "HUVEC", "K562", "NHEK", "NHLF")
urls = c("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz",
         "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmH1hescHMM.bed.gz",
         "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHepg2HMM.bed.gz",
         "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHmecHMM.bed.gz",
         "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHsmmHMM.bed.gz",
         "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHuvecHMM.bed.gz",
         "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmK562HMM.bed.gz",
         "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmNhekHMM.bed.gz",
         "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmNhlfHMM.bed.gz")

library(data.table)
ChromHMM = data.frame()
for (i in 1:length(biosamples)){
  biosample = biosamples[i]
  url = urls[i]
  
  data = data.frame(fread(url))
  set = c("4_Strong_Enhancer", "5_Strong_Enhancer", "6_Weak_Enhancer", "7_Weak_Enhancer")
  data = data[which(data$V4 %in% set), 1:4]
  data = cbind(data, biosample)
  ChromHMM = rbind(ChromHMM, data)
}
rm(data, biosample, biosamples, i, set, url, urls)

ChromHMM = ChromHMM[!duplicated(ChromHMM),]
source("./scripts/bed2coord.R")
coord = bed2coord(ChromHMM[1:3])
ChromHMM = cbind(enh_ID = paste0("chr19_", gsub("[:-]", "_", coord)), ChromHMM[1:3], coord, orig_assembly = "hg19", 
                 original_ID = NA, crossref = NA, enh_PMID = 21441907, ChromHMM[5], enh_method = "ChIP-seq,Hidden Markov Model",
                 ChromHMM[4], source = "ChromHMM")
colnames(ChromHMM)[c(2:4,12)] = c("orig_chr", "orig_start", "orig_end", "type")
rm(coord)


#2. Enhancers
enh19 = ChromHMM[1:6]
enh19 = enh19[!duplicated(enh19),]
library(parallel)
source("./scripts/hgtohgparallel.R")
hg19tohg38 = hgtohgparallel(enh19[2:5], 19, 38, detectCores()/2)

source("./scripts/size_comparison_parallel.R")
size_hg2hg = size_comparison_parallel(hg19tohg38, 0.05)

source("./scripts/filtAfterSizeComp_parallel.R")
filtered = filtAfterSizeComp_parallel(hg19tohg38, size_hg2hg, enh19[2:5])

enhancers = cbind(filtered$merge_table, score = size_hg2hg$score)
indexes = which(is.na(enhancers$new_coord))
enhancers[indexes, "score"] = NA
rm(indexes)


#3. Enrich original set
set = cbind(enhancers[4:7], current_assembly = "hg38", score = enhancers$score)
colnames(set) = c("coord", "current_chr", "current_start", "current_end", "current_assembly", "score")
ChromHMM = merge(ChromHMM, set, by = "coord")
ChromHMM = cbind(ChromHMM[c(2:6,14:17)], minimum_ratio = 0.95, ChromHMM[c(18,7:13)], hgnc_symbol_target_genes = NA,
                 enh2gene_PMID = NA, enh2gene_method = NA, hgnc_symbol_TFs = NA, TFs2enh_PMID = NA, TFs2enh_method = NA,
                 disease = NA, disease_PMID = NA, disease_method = NA, refsnp_ID = NA, mutation_PMID = NA, mutation_method = NA)
colnames(ChromHMM)[15] = "biosample_name"
set = ChromHMM
set$enh_method = "Hidden Markov Model"
ChromHMM$enh_method = "ChIP-seq"
ChromHMM = rbind(ChromHMM, set)
ChromHMM = ChromHMM[!duplicated(ChromHMM),]
rm(set, enh19)


#4. Write and save
indexes = which(is.na(ChromHMM$current_chr))
filt_ChromHMM = ChromHMM[-indexes,]
rm(indexes)

dir.create("./ChromHMM_results")
save.image("./ChromHMM_results/data_ChromHMM.RData")
write.table(filt_ChromHMM, "./ChromHMM_results/ChromHMM.tsv", col.names = T, row.names = F, quote = F, sep = "\t")