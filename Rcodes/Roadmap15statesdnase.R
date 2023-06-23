#1.Read files
#https://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_enh/
library(data.table)
path = "./Roadmap15statesdnase_data/dnase_bed_files_enh/"
list_files = list.files(path)
data = data.frame()
for (i in 1:length(list_files)){
  sample = unlist(strsplit(list_files[i], "[_|\\.]"))[3]
  file = fread(paste0(path, list_files[i]))
  file = cbind(file[,1:4], sample)
  data = rbind(data, file)
}
rm(file, i, list_files, path, sample)
data = data[!duplicated(data),]


#2. Biosamples
biosamples = fread("https://egg2.wustl.edu/roadmap/data/byFileType/metadata/EID_metadata.tab")
biosamples = biosamples[,c(1,6)]
colnames(biosamples) = c("sample", "biosample_name")
data = merge(data, biosamples, by = "sample", all.x = T)
data = data[,-1]


#3. Enhancer table
enh_ID = paste0("hg19_", gsub("[:-]", "_", data$V4))
enh_Roadmap = data.frame(enh_ID, orig_chr = data$V1, orig_start = data$V2, orig_end = data$V3, orig_assembly = "hg19",
                         minimum_ratio = "-", score = "-", original_ID = "-", crossref = "-", enh_PMID = "25693563", 
                         biosample_name = data$biosample_name, enh_method = "ChIP-seq,ChromHMM,DNase-seq", type = "-", 
                         source = "Roadmap", hgnc_symbol_target_genes = "-", enh2gene_PMID = "-", enh2gene_method = "-", 
                         hgnc_symbol_TFs = "-", TFs2enh_PMID = "-", TFs2enh_method = "-", disease = "-", disease_PMID = "-", 
                         disease_method = "-", refseq_ID = "-", mutation_PMID  = "-", mutation_method = "-")

#4. Enh19 to enh38
enh19 = data[,-5]
rm(data)
colnames(enh19) = c("orig_chr", "orig_start", "orig_end", "coord")
enh19 = cbind(enh_ID, enh19, orig_assembly = "hg19")
rm(enh_ID, biosamples)
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

enh_ID = paste0("hg19_", gsub("[:-]", "_", enhancers$orig_coord))
enhancers = cbind(enh_ID, enhancers[,1:3], orig_assembly = "hg19", enhancers[,5:7], current_assembly = "hg38", 
                  minimum_ratio = 0.95, score = enhancers$score)
indexes = which(is.na(enhancers$new_chr))
enhancers$minimum_ratio[indexes] = enhancers$current_assembly[indexes] = NA
enhancers = data.frame(enhancers, stringsAsFactors = F)
enhancers$new_chr = as.character(enhancers$new_chr)
enhancers[is.na(enhancers)] = "-"
colnames(enhancers) = c("enh_ID", "orig_chr", "orig_start", "orig_end", "orig_assembly",
                        "current_chr", "current_start", "current_end", "current_assembly", 
                        "minimum_ratio", "score")
rm(enh_ID, indexes)
enh_Roadmap = merge(enhancers[,-c(2:5)], enh_Roadmap[,-c(6:7)], by = "enh_ID", all.x = T)
enh_Roadmap$enh_method = "ChIP-seq"
set1 = set2 = enh_Roadmap
set1$enh_method = "ChromHMM"
set2$enh_method = "DNase-seq"
enh_Roadmap = rbind(enh_Roadmap, set1, set2)
rm(set1, set2)


#5. Filtering data
indexes = which(enh_Roadmap$current_chr == "-")
filt_Roadmap = enh_Roadmap[-indexes,]


#6. Saving files
dir.create("./Roadmap15statesdnase_results/")
save.image("./Roadmap15statesdnase_results/Roadmap_data.RData")
write.table(filt_Roadmap, "./Roadmap15statesdnase_results/Roadmap.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
