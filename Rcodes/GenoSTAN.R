#1. Reading files
library(data.table)
orig_GenoSTAN = fread("https://www.cmm.in.tum.de/public/paper/GenoSTAN/GenoSTAN_enhancers.bed.gz")
biosamples = fread("https://www.cmm.in.tum.de/public/paper/GenoSTAN/EID_annotation.tab")


#2. Processing table
#Merge tables
enh_GenoSTAN = orig_GenoSTAN[,1:4]
enh_GenoSTAN$V4 = unlist(lapply(enh_GenoSTAN$V4, function(biosample){unlist(strsplit(biosample, "_"))[1]}))
colnames(enh_GenoSTAN)[4] = "EID"
enh_GenoSTAN = merge(enh_GenoSTAN, biosamples[,5:6], all.x = T)
enh_GenoSTAN = enh_GenoSTAN[,-1]

#Coord and Enh_ID
source("./scripts/bed2coord.R")
coord = bed2coord(enh_GenoSTAN[,1:3])
enh_ID = paste0("hg38_", gsub("[:-]", "_", coord))
enh_GenoSTAN = data.frame(enh_ID, orig_chr = enh_GenoSTAN$V1, orig_start = enh_GenoSTAN$V2, orig_end = enh_GenoSTAN$V3, 
                          orig_assembly = "hg38", current_chr = enh_GenoSTAN$V1, current_start = enh_GenoSTAN$V2, 
                          current_end = enh_GenoSTAN$V3, current_assembly = "hg38", minimum_ratio = "-", score = "-", 
                          original_ID = "-", crossref = "-", enh_PMID = 26206277, biosample_name =  enh_GenoSTAN$`Epigenome name`, 
                          enh_method = "ChIP-seq,GenoSTAN", type = "-", source = "GenoSTAN", hgnc_symbol_target_genes = "-",
                          enh2gene_PMID = "-", enh2gene_method = "-", hgnc_symbol_TFs = "-", TFs2enh_PMID = "-", TFs2enh_method = "-",
                          disease = "-", disease_PMID = "-", disease_method = "-", refseq_ID = "-", mutation_PMID = "-", 
                          mutation_method = "-")
enh_GenoSTAN = enh_GenoSTAN[!duplicated(enh_GenoSTAN),]
enh_GenoSTAN$enh_method = "ChIP-seq"
set = enh_GenoSTAN
set$enh_method = "GenoSTAN"
enh_GenoSTAN = rbind(enh_GenoSTAN, set)
rm(enh_ID, coord, set)
enh_GenoSTAN = enh_GenoSTAN[order(enh_GenoSTAN$enh_ID),]


#3. Saving files
dir.create("./GenoSTAN_results/")
save.image("./GenoSTAN_results/GenoSTAN_data.RData")
write.table(enh_GenoSTAN, "./GenoSTAN_results/GenoSTAN.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
