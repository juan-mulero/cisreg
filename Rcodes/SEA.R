#1. Reading files
library(data.table)
orig_SEA = fread("http://218.8.241.248:8080/SEA3/download/SEA00101.bed") #Data version 2018-06-12
SEA = orig_SEA[, -c(11:15,17)]
colnames(SEA) = c("original_ID", "orig_chr", "orig_start", "orig_end", "name", "length", "biosample_name", "mean", "median",
                  "source", "external_gene_name_target_genes", "enh_method", "type", "coding")


#2. Complete table
source("./scripts/bed2coord.R")
coord = bed2coord(SEA[,2:4])
SEA = cbind(enh_ID = paste0("hg38_", gsub("[:-]", "_", coord)), SEA[,2:4], orig_assembly = "hg38", 
              current_chr = SEA$orig_chr, current_start = SEA$orig_start, current_end = SEA$orig_end, current_assembly = "hg38",
              minimum_ratio = "-", score = "-", SEA[,1], 
              crossref = paste0("http://218.8.241.248:8080/SEA3/DetailInfoAction?species=Human&SEID=", SEA$original_ID),
              enh_PMID = 31667506, SEA[,c(7,12,13,10,11)], enh2gene_PMID = 31667506, enh2gene_method = "closest", 
              hgnc_symbol_TFs = "-", TFs2enh_PMID = "-", TFs2enh_method = "-", disease = "-", disease_PMID = "-", disease_method = "-",
              refseq_ID = "-", mutation_PMID = "-", mutation_method = "-", SEA[,c(5,6,8,9,14)])


#3. Genes
source("./scripts/UpdateGeneSymbol.R")
genes = unique(SEA$external_gene_name_target_genes)
t = proc.time()
upd_genes = UpdateGeneSymbol(genes)
t = proc.time() - t

#some genes are incorrectly annotated as dates, so we correct these annotations
old_genes = c("1-Dec", "1-Mar", "10-Mar", "10-Sep", "11-Mar", "11-Sep", "12-Sep", "14-Sep", "2-Mar", "2-Sep", "3-Mar", "3-Sep",
              "4-Mar", "4-Sep", "5-Mar", "5-Sep", "6-Mar", "6-Sep", "7-Mar", "7-Sep", "8-Mar", "8-Sep", "9-Mar", "9-Sep")
new_genes = c("DELEC1","MARCHF1","MARCHF10","SEPTIN10","MARCHF11","SEPTIN11","SEPTIN12","MARCHF2","SEPTIN2","MARCHF3","SEPTIN3",
              "MARCHF4","SEPTIN4","MARCHF5","SEPTIN5","MARCHF6","SEPTIN6","MARCHF7","SEPTIN7","MARCHF8","SEPTIN8","MARCHF9","SEPTIN9")
for (i in 1:length(old_genes)){
  index = which(upd_genes$external_gene_name == old_genes[i])
  upd_genes$hgnc_symbol[index] = new_genes[i]
}
indexes = which(is.na(upd_genes$hgnc_symbol))
upd_genes$hgnc_symbol[indexes] = "-"
colnames(upd_genes)[1] = "external_gene_name_target_genes"
SEA = merge(SEA, upd_genes, by = "external_gene_name_target_genes", all.x = T)
rm(coord, genes, i, index, indexes, new_genes, old_genes)

SEA = SEA[,c(2:19,1,36,20:35)]
colnames(SEA)[20] = "hgnc_symbol_target_genes"


#4. Filtered data
filt_SEA = SEA[, -c(19,32:36)]
indexes = which(filt_SEA$hgnc_symbol_target_genes == "-")
filt_SEA$enh2gene_method[indexes] = filt_SEA$enh2gene_PMID[indexes] = "-"
rm(indexes)

chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
        "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
indexes = which(filt_SEA$orig_chr %in% chr)
filt_SEA = filt_SEA[indexes,]
rm(indexes)


#5. Save data
dir.create("./SEA_results/")
save.image("./SEA_results/SEA_data.RData")
write.table(filt_SEA, "./SEA_results/SEA.tsv", col.names = T, row.names = F, quote = F, sep = "\t")