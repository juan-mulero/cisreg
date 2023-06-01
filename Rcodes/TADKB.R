#Reading and preprocess:
library(data.table)
dir.create("./TADKB_hg19TAD")
temp = tempfile()
download.file("http://dna.cs.miami.edu/TADKB/download/TAD_annotations.tar.gz", temp)
untar(temp, exdir = "./TADKB_hg19TAD")
rm(temp)
orig_TAD_files = list.files("./TADKB_hg19TAD/TAD_annotations/TADs/")
orig_TAD = data.frame()
human_cells = c("GM12878", "HMEC", "NHEK", "IMR90", "KBM7", "K562", "HUVEC")
for (i in 1:length(orig_TAD_files)){
  split = unlist(strsplit(orig_TAD_files[i], "[_\\.]"))
  if (split[2] %in% human_cells){
    file = fread(paste0("./TADKB_hg19TAD/TAD_annotations/TADs/", orig_TAD_files[i]))
    chr = unique(file$V1)
    for (j in 1:length(chr)){
      subset = file[file$V1 == chr[j],]
      subset = data.frame(cbind(ID = seq(1,nrow(subset)), subset, biosample_name = split[2], 
                     TAD_method = split[1], TAD_submethod = split[3], resolution = split[4]))
      #Only TAD detected by Hi-C have been implemented on the website:
      if (split[1] == "HiC"){
        crossref = paste0("http://dna.cs.miami.edu/TADKB/domain.php?sp=hum&cl=", split[2], "&rg=hg19&chr=", gsub("chr", "", chr[j]), 
                        "&se=", paste(subset$V2, subset$V3, sep = "_"), "&id=", subset$ID, "&res=", split[4], "&caller=", split[3])
      } else {crossref = NA}
      subset = data.frame(cbind(subset[2:4], crossref, subset[5:8]))
      orig_TAD = rbind(orig_TAD, subset)
    }
  }
}
rm(file, subset, split, orig_TAD_files, human_cells, chr, crossref, i, j)
orig_TAD = data.frame(orig_TAD[!duplicated(orig_TAD),])
#DI --> Directionality Index
#IS --> Insulation Score (IS)
#GMAP --> Gaussian Mixture model And Proportion test

#hg19 to hg38
source("./scripts/bed2coord.R")
coord = bed2coord(orig_TAD[1:3])
assembly = "hg19"
TAD_ID = paste(assembly, gsub("[:-]", "_", coord), sep = "_")
orig_TAD = cbind(TAD_ID, orig_TAD[1:3], coord, orig_assembly = assembly, orig_TAD[4:8])
colnames(orig_TAD)[2:4] = c("orig_chr", "orig_start", "orig_end")
rm(assembly, coord, TAD_ID)

source("./scripts/hgtohg_onestep.R")
TAD19 = orig_TAD[1:6]
TAD19 = TAD19[!duplicated(TAD19),]
TAD19to38 = hgtohg_onestep(TAD19[2:5], 19, 38, 0.05)

TAD38 = TAD19to38$filtered$new_bedTable
colnames(TAD38)[4] = "coord"
TAD38 = cbind(TAD38, current_assembly = "hg38")

TAD_coord = TAD19to38$filtered$IDs
colnames(TAD_coord) = c("coord", "hg19")
TAD_coord = cbind(TAD_coord, minimum_ratio = 0.95, score = TAD19to38$size_comparison$score)

TAD38 = merge(TAD38, TAD_coord, by = "coord", all.x = T)
rm(TAD_coord)
TAD38 = TAD38[!duplicated(TAD38),]
colnames(TAD38)[c(1:4,6)] = c("hg38", "current_chr", "current_start", "current_end", "coord")
TAD38 = merge(TAD38, TAD19[c(1,5)], by = "coord", all.x = T)
TAD38 = TAD38[!duplicated(TAD38),]
TAD38 = TAD38[c(9,3:5,2,6:8)]
colnames(TAD38)[5] = "coord"

#Enrich table
TAD = merge(TAD19[c(1:4,6)], TAD38[-5], by = "TAD_ID", all = T)
TAD = TAD[!duplicated(TAD),]
TAD = merge(TAD, orig_TAD[c(1,7:11)], by = "TAD_ID", all = T)
TAD = data.frame(cbind(TAD[1:11], original_ID = NA, crossref = TAD$crossref, TAD_PMID = 30871473, TAD[13:16], type = "TAD",
                        source = "TADKB"))
set1 = TAD[-17]
set2 = TAD[-16]
colnames(set2)[16] = "TAD_method"
final_TAD = rbind(set1, set2)
final_TAD = final_TAD[order(final_TAD$orig_chr, final_TAD$orig_start, final_TAD$orig_end),]
indexes = which(is.na(final_TAD$current_chr))
if (length(indexes > 0)){filt_TAD = final_TAD[-indexes,]}
filt_TAD[is.na(filt_TAD)] = "-"

#Save files
rm(indexes)
dir.create("./TADKB_results")
save.image("./TADKB_results/data_TADKB.RData")
write.table(filt_TAD, "./TADKB_results/TADKB.tsv", quote = F, sep = "\t", row.names = F, col.names = T)
    