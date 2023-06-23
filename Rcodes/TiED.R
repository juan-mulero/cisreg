#1. Read files
library(readxl)
orig_enh_TiED = read_xlsx("./TiED_data/Active_enhancers_of_10_tissues.xlsx", skip = 2)
orig_enh_SNP = read_xlsx("./TiED_data/SNPs close to enhancers in 10 tissues.xlsx")


#2. Enhancer table
source("./scripts/bed2coord.R")
coord = bed2coord(orig_enh_TiED[,1:3])
enh_ID = paste0("hg19_", gsub("[:-]", "_", coord))
enh19 = cbind(enh_ID, orig_enh_TiED[1:3], coord, orig_assembly = "hg19")
colnames(enh19)[2:4] = c("orig_chr", "orig_start", "orig_end")
enh19 = enh19[!duplicated(enh19),]

enh_SNP = orig_enh_SNP[,c(4,9)]
colnames(enh_SNP) = c("Enhancer", "refseq_ID")

crossref = paste0("http://lcbb.swjtu.edu.cn/TiED/browseEnh/?enhancer=", orig_enh_TiED$Enhancer)
enh_TiED = cbind(orig_enh_TiED, enh_ID, crossref)
enh_TiED = merge(enh_TiED, enh_SNP, by = "Enhancer", all.x = T)
enh_TiED = data.frame(enh_ID = enh_TiED$enh_ID, orig_chr = enh_TiED$Chr, orig_start = enh_TiED$`Chr start`, 
                      orig_end = enh_TiED$`Chr end`, orig_assembly = "hg19", original_ID = enh_TiED$Enhancer, 
                      crossref = enh_TiED$crossref, enh_PMID = 30123079, biosample_name = enh_TiED$Tissues, 
                      enh_method = "ChIP-seq data of H3K27ac, H3K4me1, H3K4me3 and DHS", type = "active enhancers", 
                      source = "TiED", hgnc_symbol_target_genes = "-", enh2gene_PMID = "-", enh2gene_method = "-", 
                      hgnc_symbol_TFs = "-", TFs2enh_PMID = "-", TFs2enh_method = "-", disease = "-", disease_PMID = "-", 
                      disease_method = "-", refseq_ID = enh_TiED$refseq_ID, mutation_PMID = 30123079, mutation_method = "overlap")
enh_TiED[is.na(enh_TiED)] = "-"
indexes = which(enh_TiED$refseq_ID == "-")
enh_TiED$mutation_PMID[indexes] = enh_TiED$mutation_method[indexes] = "-"

source("./scripts/SplitColumn.R")
enh_TiED = SplitColumn(enh_TiED, 9, ",")


#3. hg19 to hg38
source("./scripts/hgtohg_onestep.R")
hg19tohg38 = hgtohg_onestep(enh19[2:6], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
hg19tohg38_tables = hgtohg2tables("hg19", "hg38", enh19, hg19tohg38)

enh_TiED = merge(hg19tohg38_tables$enhancers, enh_TiED[,c(1,6:24)], by = "enh_ID", all.x = T)
enh_TiED[is.na(enh_TiED)] = "-"


#4. Filtering and saving files
#Methods
set1 = set2 = set3 = set4 = enh_TiED
set1$enh_method = "H3K27ac"
set2$enh_method = "H3K4me1"
set3$enh_method = "H3K4me3"
set4$enh_method = "DHS"
enh_TiED = rbind(set1, set2, set3, set4)
rm(set1, set2, set3, set4, coord, crossref, enh_ID)

#Filtering
indexes = which(enh_TiED$current_chr == "-")
filt_TiED = enh_TiED[-indexes,]
filt_TiED = filt_TiED[order(filt_TiED$enh_ID),]
rm(indexes)

#Saving
dir.create("./TiED_results/")
save.image("./TiED_results/TiED_data.RData")
write.table(filt_TiED, "./TiED_results/TiED.tsv", col.names = T, row.names = F, quote = F, sep = "\t")