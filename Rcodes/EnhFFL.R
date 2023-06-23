#1. Reading files
library(data.table)
library(readxl)
enh_miRNA_gene = fread("http://lcbb.swjtu.edu.cn/EnhFFL/static_/files/Enhancer_miRNA_Gene.txt.gz")
TF_enh = fread("http://lcbb.swjtu.edu.cn/EnhFFL/static_/files/TF_Enhancer.txt.gz")
TF_enh_gene = fread("http://lcbb.swjtu.edu.cn/EnhFFL/static_/files/TF_Enhancer_Gene.txt.gz")
TF_enh_miRNA = fread("http://lcbb.swjtu.edu.cn/EnhFFL/static_/files/TF_Enhancer_miRNA.txt.gz")
colnames(enh_miRNA_gene)[1] = colnames(TF_enh)[1] = colnames(TF_enh_gene)[1] = colnames(TF_enh_miRNA)[1] = "species"

#Filter
enh_miRNA_gene = enh_miRNA_gene[enh_miRNA_gene$species == "Human",]
TF_enh = TF_enh[TF_enh$species == "Human",]
TF_enh_gene = TF_enh_gene[TF_enh_gene$species == "Human",]
TF_enh_miRNA = TF_enh_miRNA[TF_enh_miRNA$species == "Human",]


#2. Set of enhancers
#Enrich files
source("./scripts/bed2coord.R")
coord = bed2coord(enh_miRNA_gene[,c(3:5)])
enh_miRNA_gene = cbind(enh_miRNA_gene[,c(1:2)], enh_ID = paste0("hg19_", gsub("[:-]", "_", coord)), coord, enh_miRNA_gene[,3:24])
colnames(enh_miRNA_gene)[3:7] = c("enh_ID", "coord", "chr", "start", "end")

TF_enh = cbind(TF_enh[,1:7], enh_ID = paste0("hg19_", gsub("[:-]", "_", TF_enh$enh_position)), TF_enh[,8:11])
source("./scripts/coord2bed.R")
bed = coord2bed(TF_enh$enh_position)
TF_enh = cbind(TF_enh[,c(1:8)], bed, TF_enh[,c(10:12)])

bed = coord2bed(TF_enh_gene$enh_position)
TF_enh_gene = cbind(TF_enh_gene[,1:7], enh_ID = paste0("hg19_", gsub("[:-]", "_", TF_enh_gene$enh_position)), bed, TF_enh_gene[,9:15])

bed = coord2bed(TF_enh_miRNA$enh_position)
TF_enh_miRNA = cbind(TF_enh_miRNA[,1:7], enh_ID = paste0("hg19_", gsub("[:-]", "_", TF_enh_miRNA$enh_position)), bed, TF_enh_miRNA[,9:15])
rm(bed, coord)

enh19 = rbind(enh_miRNA_gene[,c(3,5:7,4)], TF_enh[,8:12], TF_enh_gene[,8:12], TF_enh_miRNA[,8:12])
enh19 = enh19[!duplicated(enh19),]
colnames(enh19)[2:4] = c("orig_chr", "orig_start", "orig_end")
enh19 = cbind(enh19, orig_assembly = "hg19")


#hg19 to hg38
source("./scripts/hgtohg_onestep.R")
enh19to38 = hgtohg_onestep(enh19[,2:5], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
tables_EnhFFL = hgtohg2tables("hg19", "hg38", enh19, enh19to38)
enh_EnhFFL = tables_EnhFFL$enhancers


#Enhancer set
set = rbind(enh_miRNA_gene[,c(3,8,2,9)], TF_enh[,c(8,13,2,14)], TF_enh_gene[,c(8,13,2,14)], TF_enh_miRNA[,c(8,13,2,14)])
set = set[!duplicated(set),]
enh_EnhFFL = merge(enh_EnhFFL, set, by = "enh_ID", all.x = T)
colnames(enh_EnhFFL)[12:14] = c("original_ID", "biosample_name", "type")
enh_EnhFFL = cbind(enh_EnhFFL[,1:12], 
                   crossref = paste0("http://lcbb.swjtu.edu.cn/EnhFFL/details/?term=", enh_EnhFFL$original_ID, "&subtype=enhancer&species=human"),
                   enh_PMID = 35694152, enh_EnhFFL[,13], enh_method = "ChIP-seq", enh_EnhFFL[,14], source = "EnhFFL")
rm(set)


#3. Genes
genes = c(enh_miRNA_gene$gene_name, TF_enh_gene$gene_name)
genes = unique(genes)
source("./scripts/UpdateGeneSymbol.R")
t = proc.time()
upd_genes = UpdateGeneSymbol(genes)
t = proc.time() - t; t #9145.53
colnames(upd_genes) = c("external_target_gene_name", "hgnc_symbol_target_genes")


#4. TF
TFs = c(TF_enh$tf_name, TF_enh_gene$tf_name, TF_enh_miRNA$tf_name)
TFs = unique(TFs)
t = proc.time()
upd_TFs = UpdateGeneSymbol(TFs)
t = proc.time() - t; t #151.53   
colnames(upd_TFs) = c("external_tf_name", "hgnc_symbol_TFs")


#5. miRNA
miRNA = c(enh_miRNA_gene$mirna_name, TF_enh_miRNA$mirna_name)
miRNA = unique(miRNA)
t = proc.time()
upd_miRNA = UpdateGeneSymbol(miRNA)
t = proc.time() - t; t #3008.76 
colnames(upd_miRNA) = c("external_miRNA_name", "hgnc_symbol_miRNA")


#6. Enrich table
enrich_table = data.frame()

#Genes:
set = enh_miRNA_gene[,c(2,3,14,23)]
colnames(set) = c("biosample_name", "enh_ID", "external_miRNA_name", "external_target_gene_name")
set = merge(set, upd_genes, by = "external_target_gene_name", all.x = T)
set = set[!duplicated(set),]
set = merge(set, upd_miRNA, by = "external_miRNA_name", all.x = T)
set = set[!duplicated(set),]
set = cbind(set[,c(4,3,2,5)], external_tf_name = NA, hgnc_symbol_TFs = NA, set[,c(1,6)])
enrich_table = rbind(enrich_table, set)

set = TF_enh_gene[,c(2,4,8,17)]
colnames(set) = c("biosample_name", "external_tf_name", "enh_ID", "external_target_gene_name")
set = merge(set, upd_genes, by = "external_target_gene_name", all.x = T)
set = set[!duplicated(set),]
set = merge(set, upd_TFs, by = "external_tf_name", all.x = T)
set = set[!duplicated(set),]
set = cbind(set[,c(4,3,2,5,1,6)], external_miRNA_name = NA, hgnc_symbol_miRNA = NA)
enrich_table = rbind(enrich_table, set)

#TFs
set = TF_enh[,c(2,4,8)]
colnames(set) = c("biosample_name", "external_tf_name", "enh_ID")
set = merge(set, upd_TFs, by = "external_tf_name", all.x = T)
set = set[!duplicated(set),]
set = cbind(set[,c(3,2)], external_target_gene_name = NA, hgnc_symbol_target_genes = NA, set[,c(1,4)], 
            external_miRNA_name = NA, hgnc_symbol_miRNA = NA)
enrich_table = rbind(enrich_table, set)

set = TF_enh_miRNA[,c(2,4,8,17)]
colnames(set) = c("biosample_name", "external_tf_name", "enh_ID", "external_miRNA_name")
set = merge(set, upd_TFs, by = "external_tf_name", all.x = T)
set = set[!duplicated(set),]
set = merge(set, upd_miRNA, by = "external_miRNA_name", all.x = T)
set = set[!duplicated(set),]
set = cbind(set[,c(4,3)], external_target_gene_name = NA, hgnc_symbol_target_genes = NA, set[,c(2,5,1,6)])
enrich_table = rbind(enrich_table, set)
rm(set)

#Complete set
enh_EnhFFL = merge(enh_EnhFFL, enrich_table, by = c("enh_ID", "biosample_name"), all.x  = T)
enh_EnhFFL = enh_EnhFFL[!duplicated(enh_EnhFFL),]
enh_EnhFFL = cbind(enh_EnhFFL[,c(1,3:15,2,16:20)], enh2gene_PMID = NA, enh2gene_method = NA, 
                   enh_EnhFFL[,c(21,22)], TFs2enh_PMID = NA, TFs2enh_method = NA,
                   disease = NA, disease_PMID = NA, disease_method = NA, refseq_ID = NA, mutation_PMID = NA, mutation_method = NA,
                   enh_EnhFFL[,c(23,24)], miRNA2enh_PMID = NA, miRNA2enh_method = NA)

indexes = which(!is.na(enh_EnhFFL$external_target_gene_name))
enh_EnhFFL$enh2gene_PMID[indexes] = 35694152
enh_EnhFFL$enh2gene_method[indexes] = "distance 100kb"

indexes = which(!is.na(enh_EnhFFL$external_tf_name))
enh_EnhFFL$TFs2enh_PMID[indexes] = 35694152
enh_EnhFFL$TFs2enh_method[indexes] = "ChIP-seq"
set = enh_EnhFFL[indexes,]
set$TFs2enh_method = "overlap"
enh_EnhFFL = rbind(enh_EnhFFL, set)
rm(set)

indexes = which(!is.na(enh_EnhFFL$external_miRNA_name))
enh_EnhFFL$miRNA2enh_PMID[indexes] = 35694152
enh_EnhFFL$miRNA2enh_method[indexes] = "distance +10kb and -1kb"
enh_EnhFFL = enh_EnhFFL[order(enh_EnhFFL$enh_ID),]


#7. Filtered data
indexes = which(enh_EnhFFL$current_chr == "-")
filt_EnhFFL = enh_EnhFFL[-indexes,c(1:18,20:22,24:32,34:36)]
indexes = which(!is.na(filt_EnhFFL$hgnc_symbol_miRNA))
set = filt_EnhFFL[indexes,c(1:18,31:33,22:30)]
colnames(set)[19:21] = c("hgnc_symbol_target_genes", "enh2gene_PMID", "enh2gene_method")
filt_EnhFFL = rbind(filt_EnhFFL[,1:30], set)

indexes = which(is.na(filt_EnhFFL$hgnc_symbol_target_genes))
filt_EnhFFL$enh2gene_PMID[indexes] = NA
filt_EnhFFL$enh2gene_method[indexes] = NA

indexes = which(is.na(filt_EnhFFL$hgnc_symbol_TFs))
filt_EnhFFL$TFs2enh_PMID[indexes] = NA
filt_EnhFFL$TFs2enh_method[indexes] = NA

filt_EnhFFL[is.na(filt_EnhFFL)] = "-"
filt_EnhFFL = filt_EnhFFL[!duplicated(filt_EnhFFL),]
rm(set, indexes)


#8.Save data
dir.create("./EnhFFL_results/")
save.image("./EnhFFL_results/EnhFFL_data.RData")
write.table(filt_EnhFFL, "./EnhFFL_results/EnhFFL.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
