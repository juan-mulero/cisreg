#1. Reading data
#Enhancers and their coordinates
library(data.table)
library(rvest)
library(xml2)
#Release 104 --> version of 8 January 2021
enhancers = fread("http://ftp.ensembl.org/pub/release-104/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz")
enhancers = data.frame(enhancers)
#unique(enhancers$V3)
#[1] "open_chromatin_region"    "TF_binding_site"          "CTCF_binding_site"        "enhancer"                
#[5] "promoter"                 "promoter_flanking_region"
enhancers = enhancers[enhancers$V3 == "enhancer",]
enhancers = enhancers[c(1,4,5,9)]
enhancers$V9 = unlist(lapply(enhancers$V9, function(enh){
  enh = unlist(strsplit(enh, ";"))[1]
  enh = unlist(strsplit(enh, ":"))[2]
}))


#Activity of the enhancers in the different epigenomes
url = "http://ftp.ensembl.org/pub/release-104/regulation/homo_sapiens/RegulatoryFeatureActivity/"
read = read_html(url)
epigenome_names = html_text(html_node(read, "pre"))
epigenome_names = unlist(strsplit(epigenome_names, "\r\n"))
epigenome_names = epigenome_names[-1]
epigenome_names = unlist(lapply(epigenome_names, function(name){unlist(strsplit(name, "/"))[1]}))

library(data.table)
activity_enh_ensembl = matrix(NA, nrow = nrow(enhancers), ncol = length(epigenome_names))
colnames(activity_enh_ensembl) = epigenome_names
rownames(activity_enh_ensembl) = enhancers$V9

t = proc.time()
for (i in 1:length(epigenome_names)){
  #Release 104 --> 8 January 2021
  file = paste0("http://ftp.ensembl.org/pub/release-104/regulation/homo_sapiens/RegulatoryFeatureActivity/", 
                epigenome_names[i], "/homo_sapiens.GRCh38.", epigenome_names[i], 
                ".Regulatory_Build.regulatory_activity.20210107.gff.gz")
  file = data.frame(fread(file))
  file = file$V9[file$V3 == "enhancer"]
  for (j in 1:length(file)){
    split = unlist(strsplit(file[j], ";"))[c(1,7)]
    split = unlist(strsplit(split, "="))[c(2,4)]
    activity_enh_ensembl[split[2], epigenome_names[i]] = split[1]
  }
}
t = proc.time() - t; t
#     user    system   elapsed 
#9569.750  2262.537 13830.733 


#2. Enhancers
enh_Ensembl = enhancers
colnames(enh_Ensembl) = c("orig_chr", "orig_start", "orig_end", "original_ID")
enh_Ensembl$orig_chr = paste0("chr", enh_Ensembl$orig_chr)
source("./scripts/bed2coord.R")
coord = bed2coord(enh_Ensembl[1:3])
crossref = paste0("https://www.ensembl.org/homo_sapiens/Regulation/Summary?rf=", enh_Ensembl$original_ID)
enh_Ensembl = cbind(enhID = paste0("hg38_", gsub("[:-]", "_", coord)), enh_Ensembl[1:3], orig_assembly = "hg38",
                    current_chr = enh_Ensembl$orig_chr, current_start = enh_Ensembl$orig_start, current_end = enh_Ensembl$orig_end,
                    current_assembly = "hg38", minimum_ratio = NA, score = NA, original_ID = enh_Ensembl$original_ID, crossref, 
                    enh_PMID = 25887522, biosample_name = NA, enh_method = "DNase-seq,ChIP-seq", type = NA, source = "Ensembl")
#enh_PMID : article of "The ensembl regulatory build"
#enh_method : "Regulatory features are inferred from publicly available experimental data sets, including:
  #Open chromatin assays (DNase-seq)
  #Histone modification assays (ChIP-seq)
  #Transcription factor binding assays (ChIP-seq)
  #Reference: Ensembl Regulatrory Build Help & Docs --> "https://www.ensembl.org/info/genome/funcgen/regulatory_features.html"


#3. Other data
enh_Ensembl = cbind(enh_Ensembl, hgnc_symbol_target_genes = NA, enh2gene_PMID = NA, enh2gene_method = NA, hgnc_symbol_TFs = NA,
                    TFs2enh_PMID = NA, TFs2enh_method = NA, disease = NA, disease_PMID = NA, disease_method = NA,
                    refsnp_ID = NA, mutation_PMID = NA, mutation_method = NA)
enh_Ensembl = enh_Ensembl[!duplicated(enh_Ensembl),]
enh_Ensembl[is.na(enh_Ensembl)] = "-"

#4. Final set and save
set1 = enh_Ensembl
set1$enh_method = "DNase-seq"
set2 = enh_Ensembl
set2$enh_method = "ChIP-seq"
split_enh_Ensembl = rbind(set1, set2)
rm(set1, set2, coord, crossref, read, epigenome_names, file, i, j, split, url)

#split_enh_Ensembl$crossref = gsub("https://www.ensembl.org/homo_sapiens/Regulation/Summary?rf=", 
#                                  "https://identifiers.org/ensembl:", split_enh_Ensembl$crossref)

dir.create("./Ensembl_results")
save.image("./Ensembl_results/Data_Ensembl.RData")
write.table(split_enh_Ensembl, "./Ensembl_results/Ensembl.tsv", quote = F, sep = "\t", col.names = T, row.names = F, 
            fileEncoding = "utf8")