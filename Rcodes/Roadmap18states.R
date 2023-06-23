#1. Read files : 
#Expanded 18-state model (6 marks, 98 epigenomes) - ChromHMM
#File: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/all_hg38lift.mnemonics.bedFiles.tgz
library(data.table)
path = "./Roadmap18states_data/all_hg38lift.mnemonics_18states.bedFiles/"
list_files = list.files(path)
data = data.frame()
for (i in 1:length(list_files)){
  biosample = unlist(strsplit(list_files[i], "_"))[1]
  file = fread(paste0(path, list_files[i]))
  file = cbind(file, biosample = biosample)
  data = rbind(data, file)
}

#States of interest: 7_EnhG1, 8_EnhG2, 9_EnhA1, 10_EnhA2, 11_EnhWk, 15_EnhBiv
data = data[data$V4 %in% c("7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2", "11_EnhWk", "15_EnhBiv"),]
data = data[!duplicated(data),]
rm(file, i, list_files, path, biosample)


#2. Biosamples
biosamples = fread("https://egg2.wustl.edu/roadmap/data/byFileType/metadata/EID_metadata.tab")
biosamples = biosamples[,c(1,6)]
colnames(biosamples) = c("biosample", "biosample_name")
data = merge(data, biosamples, by = "biosample", all.x = T)
data = data[,-1]
rm(biosamples)


#3. Enhancer table
source("./scripts/bed2coord.R")
coord = bed2coord(data[,1:3])
enh_ID = paste0("hg38_", gsub("[:-]", "_", coord))
enh_Roadmap = data.frame(enh_ID, orig_chr = data$V1, orig_start = data$V2, orig_end = data$V3, orig_assembly = "hg38",
                         current_chr = data$V1, current_start = data$V2, current_end = data$V3, current_assembly = "hg38", 
                         minimum_ratio = "-", score = "-", original_ID = "-", crossref = "-", enh_PMID = "25693563", 
                         biosample_name = data$biosample_name, enh_method = "ChIP-seq,ChromHMM", type = data$V4, 
                         source = "Roadmap", hgnc_symbol_target_genes = "-", enh2gene_PMID = "-", enh2gene_method = "-", 
                         hgnc_symbol_TFs = "-", TFs2enh_PMID = "-", TFs2enh_method = "-", disease = "-", disease_PMID = "-", 
                         disease_method = "-", refseq_ID = "-", mutation_PMID  = "-", mutation_method = "-")
rm(data, coord, enh_ID)
enh_Roadmap$enh_method = "ChIP-seq" #6 marks: H3K4me1, H3K4me3, H3K36me3, H3K27me3, H3K9me3, H3K27ac
set = enh_Roadmap
set$enh_method = "ChromHMM"
enh_Roadmap = rbind(enh_Roadmap, set)
rm(set)


#4. Saving files
dir.create("./Roadmap18states_results/")
save.image("./Roadmap18states_results/Roadmap_data.RData")
write.table(enh_Roadmap, "./Roadmap18states_results/Roadmap.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
