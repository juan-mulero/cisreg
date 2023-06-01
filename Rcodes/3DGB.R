library(data.table)
dir.create("./3DGB_hg38TAD")
temp = tempfile()
download.file("http://3dgenome.fsm.northwestern.edu/downloads/hg38.TAD.zip", temp)
unzip(temp, exdir = "./3DGB_hg38TAD")
rm(temp)
hg38TAD_files = list.files("./3DGB_hg38TAD/hg38")
hg38TAD_files = hg38TAD_files[grepl(".txt", hg38TAD_files)]
biosample_names = gsub("Bowel_Small", "BowelSmall", hg38TAD_files)
biosample_names = gsub("Cortex_DLPFC", "CortexDLPFC", biosample_names)
biosample_names = gsub("Ventricle_Right", "VentricleRight", biosample_names)
biosample_names = unlist(lapply(biosample_names, function(name){unlist(strsplit(name, "_"))[1]}))
hg38TAD = data.frame()

for (i in 1:length(hg38TAD_files)){
  data = fread(paste0("./3DGB_hg38TAD/hg38/", hg38TAD_files[i]))
  data = cbind(data, biosample_name = biosample_names[i])
  hg38TAD = rbind(hg38TAD, data)
}
hg38TAD = data.frame(hg38TAD[!duplicated(hg38TAD),])
colnames(hg38TAD)[1:3] = c("orig_chr", "orig_start", "orig_end")
source("./scripts/bed2coord.R")
coord = bed2coord(hg38TAD[1:3])
coord = gsub("[:-]", "_", coord)
assembly = "hg38"
TAD = cbind(TAD_ID = paste(assembly, coord, sep = "_"), hg38TAD[1:3], orig_assembly = assembly, current_chr = hg38TAD$orig_chr, 
             current_start = hg38TAD$orig_start, current_end = hg38TAD$orig_end, current_assembly = assembly, minimum_ratio = NA,
             score = NA, original_ID = NA, crossref = NA, TAD_PMID = 30286773, biosample_name = hg38TAD$biosample_name, 
             TAD_method = "Hi-C", resolution = NA, type = "TAD", source = "3DGB")
TAD[is.na(TAD)] = "-"
rm(assembly, biosample_names, coord, hg38TAD_files, i, data)
dir.create("./3DGB_results")
save.image("./3DGB_results/data_3DGB.RData")
write.table(TAD, "./3DGB_results/3DGB.tsv", quote = F, sep = "\t", row.names = F, col.names = T)