#1. RegElements Double Elite
#1.1. Read files
path = "./GeneHancer_data/geneHancerRegElementsDoubleElite/"
files_RegElem_DE = list.files(path)
RegElem_DE = data.frame()
for (i in 1:length(files_RegElem_DE)){
  file = read.delim(paste0(path, files_RegElem_DE[i]))
  file = file[,c(1:5,10,11)]
  RegElem_DE = rbind(RegElem_DE, file)
}
rm(file, files_RegElem_DE, i, path)

#1.2. Table
source("./scripts/bed2coord.R")
coord = bed2coord(RegElem_DE[,1:3])
enh_ID = paste0("hg38_", gsub("[:-]", "_", coord))
RegElem_DE = cbind(enh_ID, RegElem_DE[,1:3], orig_assembly = "hg38", RegElem_DE[,1:3], current_assembly = "hg38",
                   minimum_ratio = NA, score = NA, original_ID = RegElem_DE$name, crossref = NA, enh_PMID = 28605766,
                   biosample_name = NA, enh_method = "automatic evidence combination", type = RegElem_DE$elementType,
                   source = "GeneHancer", RegElem_DE[,5:6], eliteness = "Double Elite")
colnames(RegElem_DE)[c(2:4,6:8,19)] = c("orig_chr", "orig_start", "orig_end", "current_chr", "current_start", "current_end",
                                           "orig_score")


#2. Interactions Double Elite
#2.1. Read files
path = "./GeneHancer_data/geneHancerInteractionsDoubleElite/"
files_Int_DE = list.files(path)
Int_DE = data.frame()
for (i in 1:length(files_Int_DE)){
  file = read.delim(paste0(path, files_Int_DE[i]))
  file = file[,c(12,17,7,5,6)]
  Int_DE = rbind(Int_DE, file)
}
rm(file, files_Int_DE, i, path)

#2.2. Output table
colnames(Int_DE) = c("original_ID", "external_gene_name", "enh2gene_method", "enh2gene_score", "enh2gene_value")
RegElem_DE = merge(RegElem_DE, Int_DE, by = "original_ID", all = T)
rm(Int_DE)


#3. RegElements
#3.1. Read files
path = "./GeneHancer_data/geneHancerRegElements/"
files_RegElem = list.files(path)
RegElem = data.frame()
for (i in 1:length(files_RegElem)){
  file = read.delim(paste0(path, files_RegElem[i]))
  file = file[,c(1:5,10,11)]
  RegElem = rbind(RegElem, file)
}
rm(file, files_RegElem, i, path)

#3.2. Table
coord = bed2coord(RegElem[,1:3])
enh_ID = paste0("hg38_", gsub("[:-]", "_", coord))
RegElem = cbind(enh_ID, RegElem[,1:3], orig_assembly = "hg38", RegElem[,1:3], current_assembly = "hg38",
                   minimum_ratio = NA, score = NA, original_ID = RegElem$name, crossref = NA, enh_PMID = 28605766,
                   biosample_name = NA, enh_method = "automatic evidence combination", type = RegElem$elementType,
                   source = "GeneHancer", RegElem[,5:6], eliteness = "-")
colnames(RegElem)[c(2:4,6:8,19)] = c("orig_chr", "orig_start", "orig_end", "current_chr", "current_start", "current_end",
                                        "orig_score")

#4. Interactions
#4.1. Read files
path = "./GeneHancer_data/geneHancerInteractions/"
files_Int = list.files(path)
Int = data.frame()
for (i in 1:length(files_Int)){
  file = read.delim(paste0(path, files_Int[i]))
  file = file[,c(12,17,7,5,6)]
  Int = rbind(Int, file)
}
rm(file, files_Int, i, path)

#4.2. Output table
colnames(Int) = c("original_ID", "external_gene_name", "enh2gene_method", "enh2gene_score", "enh2gene_value")
RegElem = merge(RegElem, Int, by = "original_ID", all = T)
rm(Int, coord, enh_ID)


#5. Merge tables
indexes = grep(",", RegElem$evidenceSources)
RegElem$eliteness[indexes] = "Elite"
indexes = which(grepl(",", RegElem$evidenceSources) & grepl(",", RegElem$enh2gene_method))
RegElem$eliteness[indexes] = "Double Elite"
rm(indexes)
RegElem[is.na(RegElem)] = "-"


#6. Enhancers
enh_GeneHancer = RegElem[RegElem$type %in% c("Enhancer", "Promoter/Enhancer"),]
source("./scripts/SplitColumn.R")
enh_GeneHancer = SplitColumn(enh_GeneHancer, 23, ",")
enh_GeneHancer = cbind(enh_GeneHancer[,c(2:12,1,13:18,22)], enh2gene_PMID = 28605766, enh2gene_method = enh_GeneHancer$enh2gene_method, 
                       hgnc_symbol_TFs = "-", TFs2enh_PMID = "-", TFs2enh_method = "-", 
                       disease = "-", disease_PMID = "-", disease_method = "-",
                       refseq_ID = "-", mutation_PMID = "-", mutation_method = "-", enh_GeneHancer[,c(19:21,24:25)]) 

#7. Update genes
indexes = grep("^ENSG[0-9]{11}$", enh_GeneHancer$external_gene_name)
ensembl_genes = unique(enh_GeneHancer$external_gene_name[indexes])
library(biomaRt)
ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes = getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol"), filters = "ensembl_gene_id",
              values = ensembl_genes, mart = ensembl)
ensembl_genes = data.frame(ensembl_gene_id = ensembl_genes)
ensembl_genes = merge(ensembl_genes, genes, by = "ensembl_gene_id", all = T)

source("./scripts/UpdateGenes.R")
genes = unique(enh_GeneHancer$external_gene_name[-indexes])
indexes = which(genes == "-")
genes = genes[-indexes]
t = proc.time()
upd_genes = UpdateGenes(genes)
t = proc.time() - t
#user  system elapsed
#3210.30    3.38 3231.92 

genes = data.frame(external_gene_name = genes)
genes = merge(genes, upd_genes, by = "external_gene_name", all.x = T)
colnames(genes) = c("external_target_gene_name", "hgnc_symbol_target_genes")


#8. Final set
colnames(enh_GeneHancer)[19] = "external_target_gene_name"
enh_GeneHancer = merge(enh_GeneHancer, genes, by = "external_target_gene_name", all.x = T)
enh_GeneHancer = enh_GeneHancer[,c(2:19,1,36,20:35)]
enh_GeneHancer[is.na(enh_GeneHancer)] = "-"

filt_GeneHancer = enh_GeneHancer[,-c(19,32:36)]
indexes = which(filt_GeneHancer$hgnc_symbol_target_genes == "-")
filt_GeneHancer$enh2gene_method[indexes] = filt_GeneHancer$enh2gene_PMID[indexes] = "-"

crossref = c()
for (i in 1:nrow(filt_GeneHancer)){
  chr = filt_GeneHancer$orig_chr[i]
  chr = gsub("chr", "", chr)
  link = paste0("https://genecards.weizmann.ac.il/geneloc-bin/display_map.pl?chr_nr=", chr,
                "&range_type=gh_id&gh_id=", filt_GeneHancer$original_ID[i], "#", filt_GeneHancer$original_ID[i])
  crossref[i] = link
}
filt_GeneHancer$crossref = enh_GeneHancer$crossref = crossref
rm(chr, crossref, i, indexes, link, upd_genes)


#9. Save files
dir.create("./GeneHancer_results/")
save.image("./GeneHancer_results/GeneHancer_data.RData")
write.table(filt_GeneHancer, "./GeneHancer_results/GeneHancer.tsv", quote = F, col.names = T, row.names = F, sep = "\t")
