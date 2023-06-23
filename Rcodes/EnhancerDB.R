#1. Enhancers
##Reading files
library(data.table)
orig_EnhancerDB = data.frame(fread("http://lcbb.swjtu.edu.cn/EnhancerDB/_download/enhancer.gz"))

##Set of enhancers and metadata
enh_EnhancerDB = orig_EnhancerDB
colnames(enh_EnhancerDB) = c("orig_chr", "orig_start", "orig_end", "original_ID", "biosample_name", "R1", "R2", "R3", "score")
source("./scripts/bed2coord.R")
coord = bed2coord(enh_EnhancerDB[1:3])
enh_ID = paste0("hg19_", gsub("[:-]", "_", coord))
crossref = paste0("http://lcbb.swjtu.edu.cn/EnhancerDB/enhancer/details?eid=", enh_EnhancerDB$original_ID)
enh_EnhancerDB = cbind(enh_ID, enh_EnhancerDB[1:3], orig_assembly = "hg19", enh_EnhancerDB[4], crossref, enh_PMID = 30689845,
                       enh_EnhancerDB[5], enh_method = "H3K4me1", type = NA, source = "EnhancerDB")
#EnhancerDB contains 65423 enhancers from FANTOM5, but these have been incorrectly annotated as VISTA enhancers. 
which(enh_EnhancerDB$original_ID == "vista1") #461689 527112
enh_EnhancerDB$enh_PMID[461689:527111] = 24670763
enh_EnhancerDB$enh_method[461689:527111] = "CAGE"

enh_EnhancerDB$enh_PMID[527112:nrow(enh_EnhancerDB)] = 17130149
enh_EnhancerDB$enh_method[527112:nrow(enh_EnhancerDB)] = "sequence conservation"
set = enh_EnhancerDB[527112:nrow(enh_EnhancerDB),]
set$enh_method = "reporter gene assay"
enh_EnhancerDB = rbind(enh_EnhancerDB, set)

indexes = which(grepl("enh", enh_EnhancerDB$original_ID))
set2 = enh_EnhancerDB[indexes,]
set2$enh_method = "H3K4me3"
set3 = enh_EnhancerDB[indexes,]
set3$enh_method = "H3K27ac"
set4 = enh_EnhancerDB[indexes,]
set4$enh_method = "DNase-seq"
enh_EnhancerDB = rbind(enh_EnhancerDB, set2, set3, set4)
rm(set, set2, set3, set4)

##hg19 to hg38
coord = bed2coord(enh_EnhancerDB[1:3])
source("./scripts/hgtohg_onestep.R")
hg19 = cbind(enh_EnhancerDB[1:4], coord, enh_EnhancerDB[5])
hg19 = hg19[!duplicated(hg19),]
hg19tohg38 = hgtohg_onestep(hg19[2:5], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
tables_hg19tohg19 = hgtohg2tables("hg19", "hg38", hg19, hg19tohg38)
enh_EnhancerDB = merge(tables_hg19tohg19$enhancers, enh_EnhancerDB[c(1,6:12)], by = "enh_ID")

enh_EnhancerDB = cbind(enh_EnhancerDB, hgnc_symbol_target_genes = NA, enh2gene_PMID = NA, enh2gene_method = NA, hgnc_symbol_TFs = NA,
                       TFs2enh_PMID = NA, TFs2enh_method = NA, disease = NA, disease_PMID = NA, disease_method = NA, refseq_ID = NA,
                       mutation_PMID = NA, mutation_method = NA)
enh_EnhancerDB = enh_EnhancerDB[!duplicated(enh_EnhancerDB),]
enh_EnhancerDB[is.na(enh_EnhancerDB)] = "-"

indexes = which(enh_EnhancerDB$current_chr == "-")
filt_EnhancerDB = enh_EnhancerDB[-indexes,]

dir.create("./EnhancerDB_results")
save.image("./EnhancerDB_results/Data_EnhancerDB.RData")
write.table(filt_EnhancerDB, "./EnhancerDB_results/EnhancerDB.tsv", col.names = T, row.names = F, quote = F, sep = "\t")


#2. Enhancer-miRNA
miRNA = data.frame(fread("http://lcbb.swjtu.edu.cn/EnhancerDB/_download/enhancer-miRNA.gz"))
colnames(miRNA) = c("orig_chr", "orig_start", "orig_end", "original_ID", "miRNA_name", "chr", "start", "end", "TSS", "strand", "qValue")
#We will not use the original IDs because for FANTOM enhancers are incorrect.
coord = bed2coord(miRNA[1:3])
enh_ID = paste0("hg19_", gsub("[:-]", "_", coord))
enh2miRNA = cbind(enh_ID, miRNA, enh2miRNA_PMID = 30689845, enh2miRNA_method = "distance 100kb")
enh2miRNA = enh2miRNA[!duplicated(enh2miRNA),]

miRNAs = unique(enh2miRNA$miRNA_name)
source("./scripts/UpdateGeneSymbol.R")
t = proc.time()
upd_miRNAs = UpdateGeneSymbol(miRNAs)
t = proc.time() - t
colnames(upd_miRNAs)[1] = colnames(miRNA)[5] = "external_miRNA_name"
colnames(upd_miRNAs)[2] = "hgnc_symbol_miRNA"
miRNA = merge(miRNA, upd_miRNAs, by = "external_miRNA_name", all.x = T)
miRNA = miRNA[!duplicated(miRNA),]
miRNA = miRNA[,c(2:5,1,12,6:11)]

colnames(enh2miRNA)[6] = "external_miRNA_name"
enh2miRNA = merge(enh2miRNA, upd_miRNAs, by = "external_miRNA_name", all.x = T)
enh2miRNA = enh2miRNA[!duplicated(enh2miRNA),]
enh2miRNA = enh2miRNA[,c(2:6,1,15,13:14,7:12)]

save.image("./EnhancerDB_results/Data_EnhancerDB.RData")
write.table(filt_EnhancerDB, "./EnhancerDB_results/EnhancerDB.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

