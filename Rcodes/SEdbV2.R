#1. Super-Enhancers
#1.1. Read data
library(data.table)
SEdb2.0 = read.table("./SEdb2.0_data/SE_package_hg38.bed", fill = T, header = T, stringsAsFactors = F)
SEdb2.0 = SEdb2.0[,c(1:6,13:21)]


#1.2. Data processing
#Genes to unique variable (the methods don't correspond with the present in the web)
new_SEdb2.0 = data.frame()
for (i in 1:nrow(SEdb2.0)){
  row = SEdb2.0[i,1:8]
  genes = unlist(SEdb2.0[i,9:15])
  indexes = which(is.na(genes) | genes == "")
  if (length(indexes) > 0){genes = genes[-indexes]}
  genes = unlist(strsplit(genes, ","))
  genes = unique(genes)
  if (length(genes) > 0){
    for (j in 1:length(genes)){
      gene = genes[j]
      new_SEdb2.0 = rbind(new_SEdb2.0, cbind(row, gene))
    }
  } else {
    gene = NA
    new_SEdb2.0 = rbind(new_SEdb2.0, cbind(row, gene))
  }
}
orig_SEdb2.0 = SEdb2.0
SEdb2.0 = new_SEdb2.0
rm(new_SEdb2.0, row, gene, genes, i, indexes, j)

#Crossref
SEdb2.0 = cbind(SEdb2.0, crossref = paste0("http://www.licpathway.net/sedb/search/search_se.php?se_id=", SEdb2.0$se_id))

#Bed to coordinates
source("./scripts/bed2coord.R")
coord = bed2coord(SEdb2.0[,c(3:5)])
enh_ID = paste0("hg38_", gsub("[:-]", "_", coord))
SEdb2.0 = cbind(enh_ID, SEdb2.0)
rm(coord, enh_ID)

#Biosamples
biosamples = read.csv("./SEdb2.0_data/SEdb_2.0_sample_information.csv", stringsAsFactors = F)
biosamples = biosamples[biosamples$Species == "Human",]
colnames(SEdb2.0)[2] = "Sample.ID"
SEdb2.0 = merge(SEdb2.0, biosamples[,c(1,3,6)], by = "Sample.ID", all.x = T)
SEdb2.0 = SEdb2.0[!duplicated(SEdb2.0),]

#Final set
SEdb2.0 = cbind(enh_ID = SEdb2.0$enh_ID, orig_chr = SEdb2.0$se_chr, orig_start = SEdb2.0$se_start, 
                     orig_end = SEdb2.0$se_end, orig_assembly = "hg38", current_chr = SEdb2.0$se_chr, 
                     current_start = SEdb2.0$se_start, current_end = SEdb2.0$se_end, current_assembly = "hg38", 
                     minimum_ratio = NA, score = NA, original_ID = SEdb2.0$se_id, crossref = SEdb2.0$crossref,
                     enh_PMID = 30371817, biosample_name = SEdb2.0$Biosample.name, enh_method = "H3K27ac", type = "SE", 
                     source = "SEdb", external_target_gene_names = SEdb2.0$gene, enh2gene_PMID = 30371817, 
                     enh2gene_method = NA, hgnc_symbol_TFs = NA, TFs2enh_PMID = NA, TFs2enh_method = NA, 
                     disease = NA, disease_PMID = NA, disease_method = NA, refseq_ID = NA, mutation_PMID = NA, 
                     mutation_method = NA, SE_rank = SEdb2.0$se_rank, ChIP_density_case = SEdb2.0$se_conser_score, 
                     ChIP_density_control = SEdb2.0$se_cas_value, original_source = SEdb2.0$Data.source)
SEdb2.0 = as.data.frame(SEdb2.0)
SEdb2.0[is.na(SEdb2.0)] = "-"
SEdb2.0 = SEdb2.0[!duplicated(SEdb2.0),]


#1.3. Gene updating
genes = unique(SEdb2.0$external_target_gene_names)
index = which(genes == "-")
genes = genes[-index]
rm(index)
source("./scripts/UpdateGenes.R")
upd_genes = UpdateGenes(genes)

genes = data.frame(external_target_gene_names = genes)
colnames(upd_genes) = c("external_target_gene_names", "hgnc_symbol_target_genes")
genes = merge(genes, upd_genes, by = "external_target_gene_names", all.x = T)
genes = genes[!duplicated(genes),]
genes[is.na(genes)] = "-"

SEdb2.0 = merge(SEdb2.0, genes, by = "external_target_gene_names", all.x = T)
SEdb2.0 = SEdb2.0[,c(2:19,1,35,20:34)]
SEdb2.0[is.na(SEdb2.0)] = "-"


#1.4. Filtered data and save files
filt_SEdb2.0 = SEdb2.0[,-c(19,32:35)]
indexes = which(filt_SEdb2.0$hgnc_symbol_target_genes == "-")
filt_SEdb2.0$enh2gene_PMID[indexes] = "-"
standard_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                 "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
indexes = which(filt_SEdb2.0$current_chr %in% standard_chr)
filt_SEdb2.0 = filt_SEdb2.0[indexes,]
rm(indexes)

dir.create("./SEdb2.0_results/")
save.image("./SEdb2.0_results/SEdbV2_data.RData")
write.table(filt_SEdb2.0, "./SEdb2.0_results/SE_SEdb.tsv", col.names = T, row.names = F, quote = F, sep = "\t")


#2. Super enhancer elements (constituent enhancers)
#2.1. Read data
const_SEdb2.0 = read.table("./SEdb2.0_data/SE_ele_package_hg38.bed", fill = T, header = T, stringsAsFactors = F)
const_SEdb2.0 = const_SEdb2.0[,1:5]


#2.2. Processing data
coord = bed2coord(const_SEdb2.0[,c(2:4)])
enh_ID = paste0("hg38_", gsub("[:-]", "_", coord))
const_SEdb2.0 = data.frame(enh_ID, orig_chr = const_SEdb2.0$ele_chr, orig_start = const_SEdb2.0$ele_start, 
                           orig_end = const_SEdb2.0$ele_end, orig_assembly = "hg38", current_chr = const_SEdb2.0$ele_chr,
                           current_start = const_SEdb2.0$ele_start, current_end = const_SEdb2.0$ele_end, 
                           current_assembly = "hg38", minimum_ratio = "-", score = "-", original_ID = "-", crossref = "-", 
                           enh_PMID = 30371817, biosample_name = "-", enh_method = "H3K27ac", type = "constituent enhancer", 
                           source = "SEdb", hgnc_symbol_target_genes = "-", enh2gene_PMID = "-", enh2gene_method = "-", 
                           hgnc_symbol_TFs = "-", TFs2enh_PMID = "-", TFs2enh_method = "-", disease = "-", disease_PMID = "-", 
                           disease_method = "-", refseq_ID = "-", mutation_PMID = "-", mutation_method = "-", 
                           Elem_rank = const_SEdb2.0$ele_score, part_of = const_SEdb2.0$se_id)

filt_const_SEdb2.0 = const_SEdb2.0[, -c(31:32)]
indexes = which(filt_const_SEdb2.0$current_chr %in% standard_chr)
filt_const_SEdb2.0 = filt_const_SEdb2.0[indexes,]
rm(indexes, coord, enh_ID)


#2.3. Saving files
save.image("./SEdb2.0_results/SEdbV2_data.RData")
write.table(filt_const_SEdb2.0, "./SEdb2.0_results/const_SEdb.tsv", col.names = T, row.names = F, quote = F, sep = "\t")


#3. Unified file
filt_SEdb = rbind(filt_SEdb2.0, filt_const_SEdb2.0)
write.table(filt_SEdb, "./SEdb2.0_results/SEdb.tsv", col.names = T, row.names = F, quote = F, sep = "\t")