#1. Data collection
library(data.table)
ENCODE = fread("http://acgt.cs.tau.ac.il/focs/data/encode_interactions.hg38.txt")
Roadmap = fread("http://acgt.cs.tau.ac.il/focs/data/roadmap_interactions.hg38.txt")
FANTOM5 = fread("http://acgt.cs.tau.ac.il/focs/data/fantom_interactions.hg38.txt")
Groseq = fread("http://acgt.cs.tau.ac.il/focs/data/groseq_interactions.hg38.txt")


#2. Enhancer data tables
#2.1. ENCODE
colnames(ENCODE) = colnames(FANTOM5) = colnames(Groseq) = colnames(Roadmap) = 
  c("p_chr", "p_start", "p_end", "p_strand", "e_chr", "e_start", "e_end", "e_strand", "p_ID", "e_ID", "distance",
    "external_gene_name", "ncbi_gene", "ensembl_ID", "refseq_ID")
source("./scripts/bed2coord.R")
coord = bed2coord(ENCODE[,c(5:7)])
enh_ID = paste0("hg38_", gsub("[:-]", "_", coord))
enh_ENCODE = cbind(enh_ID, ENCODE[,c(5:7)], orig_assembly = "hg38", ENCODE[,c(5:7)], current_assembly = "hg38", 
                   minimum_ratio = NA, score = NA, original_ID = ENCODE$e_ID, crossref = NA, enh_PMID = 22955617,
                   biosample_name = NA, enh_method = "DNase-seq", type = NA, source = "FOCS")
colnames(enh_ENCODE)[c(2:4,6:8)] = c("orig_chr", "orig_start", "orig_end", "current_chr", "current_start", "current_end")

ENCODE = cbind(ENCODE[,c(1:4)], enh_ID, ENCODE[,c(5:12)])
source("./scripts/SplitColumn.R")
ENCODE = SplitColumn(ENCODE, 13, ";")
ENCODE[ENCODE == ""] = NA
ENCODE = ENCODE[!duplicated(ENCODE),]


#2.2. Roadmap
coord = bed2coord(Roadmap[,c(5:7)])
enh_ID = paste0("hg38_", gsub("[:-]", "_", coord))
enh_Roadmap = cbind(enh_ID, Roadmap[,c(5:7)], orig_assembly = "hg38", Roadmap[,c(5:7)], current_assembly = "hg38", 
                   minimum_ratio = NA, score = NA, original_ID = Roadmap$e_ID, crossref = NA, enh_PMID = 25693563,
                   biosample_name = NA, enh_method = "DNase-seq", type = NA, source = "FOCS")
colnames(enh_Roadmap)[c(2:4,6:8)] = c("orig_chr", "orig_start", "orig_end", "current_chr", "current_start", "current_end")
Roadmap = cbind(Roadmap[,c(1:4)], enh_ID, Roadmap[,c(5:12)])
Roadmap = SplitColumn(Roadmap, 13, ";")
Roadmap[Roadmap == ""] = NA
Roadmap = Roadmap[!duplicated(Roadmap),]

set1 = set2 = enh_Roadmap
set1$enh_method = "ChIP-seq"
set2$enh_method = "Hidden Markov Model"
enh_Roadmap = rbind(enh_Roadmap, set1, set2)
rm(set1, set2)


#2.3. FANTOM5
coord = bed2coord(FANTOM5[,c(5:7)])
enh_ID = paste0("hg38_", gsub("[:-]", "_", coord))
enh_FANTOM5 = cbind(enh_ID, FANTOM5[,c(5:7)], orig_assembly = "hg38", FANTOM5[,c(5:7)], current_assembly = "hg38", 
                    minimum_ratio = NA, score = NA, original_ID = FANTOM5$e_ID, crossref = NA, enh_PMID = 24670763,
                    biosample_name = NA, enh_method = "CAGE", type = NA, source = "FOCS")
colnames(enh_FANTOM5)[c(2:4,6:8)] = c("orig_chr", "orig_start", "orig_end", "current_chr", "current_start", "current_end")
FANTOM5 = cbind(FANTOM5[,c(1:4)], enh_ID, FANTOM5[,c(5:12)])
FANTOM5 = SplitColumn(FANTOM5, 13, ";")
FANTOM5[FANTOM5 == ""] = NA
FANTOM5 = FANTOM5[!duplicated(FANTOM5),]


#2.4. Gro-seq
coord = bed2coord(Groseq[,c(5:7)])
enh_ID = paste0("hg38_", gsub("[:-]", "_", coord))
enh_Groseq = cbind(enh_ID, Groseq[,c(5:7)], orig_assembly = "hg38", Groseq[,c(5:7)], current_assembly = "hg38", 
                    minimum_ratio = NA, score = NA, original_ID = Groseq$e_ID, crossref = NA, enh_PMID = 29716618,
                    biosample_name = NA, enh_method = "GRO-seq", type = NA, source = "FOCS")
colnames(enh_Groseq)[c(2:4,6:8)] = c("orig_chr", "orig_start", "orig_end", "current_chr", "current_start", "current_end")
Groseq = cbind(Groseq[,c(1:4)], enh_ID, Groseq[,c(5:12)])
Groseq = SplitColumn(Groseq, 13, ";")
Groseq[Groseq == ""] = NA
Groseq = Groseq[!duplicated(Groseq),]

rm(coord, enh_ID)


#3. Genes
genes = unique(c(ENCODE$external_gene_name, Roadmap$external_gene_name, FANTOM5$external_gene_name, Groseq$external_gene_name))
index = which(is.na(genes))
genes = genes[-index]
rm(index)
source("./scripts/UpdateGeneSymbol.R")
t = proc.time()
upd_genes = UpdateGeneSymbol(genes)
t = proc.time() - t; t
#user   system  elapsed 
#94.68    31.28 57873.03 


#3.1. ENCODE
ENCODE = merge(ENCODE, upd_genes, by = "external_gene_name", all.x = T)
ENCODE = ENCODE[,c(2:13,1,14)]
enh_ENCODE = merge(data.frame(enh_ENCODE), ENCODE[,c(5,13,14)], by = "enh_ID", all.x = T)
enh_ENCODE = enh_ENCODE[!duplicated(enh_ENCODE),]


#3.2. Roadmap
Roadmap = merge(Roadmap, upd_genes, by= "external_gene_name", all.x = T)
Roadmap = Roadmap[,c(2:13,1,14)]
enh_Roadmap = merge(data.frame(enh_Roadmap), Roadmap[,c(5,13,14)], by = "enh_ID", all.x = T)
enh_Roadmap = enh_Roadmap[!duplicated(enh_Roadmap),]


#3.3. FANTOM5
FANTOM5 = merge(FANTOM5, upd_genes, by= "external_gene_name", all.x = T)
FANTOM5 = FANTOM5[,c(2:13,1,14)]
enh_FANTOM5 = merge(data.frame(enh_FANTOM5), FANTOM5[,c(5,13,14)], by = "enh_ID", all.x = T)
enh_FANTOM5 = enh_FANTOM5[!duplicated(enh_FANTOM5),]


#3.4. Gro-seq
Groseq = merge(Groseq, upd_genes, by= "external_gene_name", all.x = T)
Groseq = Groseq[,c(2:13,1,14)]
enh_Groseq = merge(data.frame(enh_Groseq), Groseq[,c(5,13,14)], by = "enh_ID", all.x = T)
enh_Groseq = enh_Groseq[!duplicated(enh_Groseq),]


#4. Complete set
enh_FOCS = rbind(enh_ENCODE, enh_Roadmap, enh_FANTOM5, enh_Groseq)
colnames(enh_FOCS)[19:20] = c("external_target_gene_name", "hgnc_symbol_target_genes")
enh_FOCS = cbind(enh_FOCS, enh2gene_PMID = 29716618, enh2gene_method = "FOCS", hgnc_symbol_TFs = NA, TFs2enh_PMID = NA,
                 TFs2enh_method = NA, disease = NA, disease_PMID = NA, disease_method = NA, refseq_ID = NA, 
                 mutation_PMID = NA, mutation_method = NA)
enh_FOCS[is.na(enh_FOCS)] = "-"


#5. Filtered table
filt_FOCS = enh_FOCS[,-19]
indexes = which(filt_FOCS$hgnc_symbol_target_genes == "-")
filt_FOCS$enh2gene_PMID[indexes] = filt_FOCS$enh2gene_method[indexes] = "-"
filt_FOCS = filt_FOCS[!duplicated(filt_FOCS),]

chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
        "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
indexes = which(filt_FOCS$orig_chr %in% chr)
filt_FOCS = filt_FOCS[indexes,]
rm(indexes)


#6. Save data
dir.create("./FOCS_results/")
save.image("./FOCS_results/FOCS_data.RData")
write.table(filt_FOCS, "./FOCS_results/FOCS.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
