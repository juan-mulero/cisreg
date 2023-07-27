# 1. Reading:
library(data.table)
orig_CancerEnD = data.frame(fread("https://webs.iiitd.edu.in/raghava/cancerend/cancerend.csv"))
crossref = paste0("https://webs.iiitd.edu.in/raghava/cancerend/display.php?details=", orig_CancerEnD$ID)
CancerEnD = cbind(orig_CancerEnD, crossref)

# 2. Genes
#The dataset does not have information about chromosomes. The genes were identified through the nearest gene method, 
#so we identified the chromosomes on which the genes are located.
source("./scripts/UpdateGeneSymbol.R")
genes = unique(orig_CancerEnD$Gene)
t = proc.time()
upd_genes = UpdateGeneSymbol(genes)
t = proc.time() - t
indexes = which(is.na(upd_genes$hgnc_symbol))
upd_genes$external_gene_name[indexes]
genes_NA = unique(upd_genes$external_gene_name[indexes])
# [1] "Dec-01" "Mar-01" "Mar-08" "Mar-11" "Sep-01" "Sep-02" "Sep-03" "Sep-05" "Sep-12" "Sep-14"
#We have genes incorrectly annotated as dates, so we replace them
genes_NA = c("Dec-01", "Mar-01", "Mar-08", "Mar-11", "Sep-01", "Sep-02", "Sep-03", "Sep-05", "Sep-12", "Sep-14")
new_genes = c("DELEC1", "MARCHF1", "MARCHF8", "MARCHF11", "SEPTIN1", "SEPTIN2", "SEPTIN3", "SEPTIN5", "SEPTIN12", "SEPTIN14")
colnames(CancerEnD) = c("original_ID", "disease", "orig_start", "orig_end", "hasMutation", "Reference.allele", "Mutated.allele",
                        "COSMIC_ID", "CNV_type", "Survival", "external_gene_name_target_genes", "gene_start", "gene_end", 
                        "source", "crossref")
for (i in 1:length(genes_NA)){
  indexes = which(upd_genes$external_gene_name == genes_NA[i])
  upd_genes$external_gene_name[indexes] = upd_genes$hgnc_symbol[indexes] = new_genes[i]
  
  indexes = which(CancerEnD$external_gene_name_target_genes == genes_NA[i])
  CancerEnD$external_gene_name_target_genes[indexes] = new_genes[i]
}

colnames(upd_genes) = c("external_gene_name_target_genes", "hgnc_symbol_target_genes")
CancerEnD = merge(CancerEnD, upd_genes, by = "external_gene_name_target_genes")
CancerEnD = CancerEnD[!duplicated(CancerEnD),]
rm(genes_NA, i, indexes, new_genes, genes)


#Chromosomes
library(biomaRt)
genes = unique(upd_genes$hgnc_symbol_target_genes)
ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
chr = getBM(attributes = c("chromosome_name", "hgnc_symbol"), filters = "hgnc_symbol", values = genes, mart = ensembl)
indexes = grep("^[0-9]|^X|^Y$", chr$chromosome_name)
chr = chr[indexes,]
chr$chromosome_name = paste0("chr", chr$chromosome_name)
colnames(chr) = c("orig_chr", "hgnc_symbol_target_genes")
CancerEnD = merge(CancerEnD, chr, by = "hgnc_symbol_target_genes", all.x = T)
unique(CancerEnD$hgnc_symbol_target_genes[is.na(CancerEnD$orig_chr)]) #"AKAP2" "PALM2" "SPHAR"
CancerEnD$orig_chr[CancerEnD$hgnc_symbol_target_genes == "AKAP2"] = "chr9"
CancerEnD$orig_chr[CancerEnD$hgnc_symbol_target_genes == "PALM2"] = "chr9"
CancerEnD$orig_chr[CancerEnD$hgnc_symbol_target_genes == "SPHAR"] = "chr1"
rm(genes, indexes)

#The enhancers were taken from study PMID: 29625054 which are mapped to the hg19 reference genome. 
#In turn, these sequences belong to FANTOM5 (PMID 24670763).
#PMID: 29625054 --> "we present the detection and characterization of a large number of expressed enhancers 
#in a genome-wide analysis of 8928 tumor samples across 33 cancer types using TCGA RNA-seq data".
#Therefore, it is important to consider that diseases are the biological sample, because it is not a relationship between 
#the enhancer and the disease, but simply an expression of the enhancer in the disease (it can be normal, good or bad).
#These genes were linked to genes by distance (0.1 MB window) and subsequently by eQTL (PMID 32360910, from CancerEnD)
#However, the coordinates of the genes are annotated with the hg38 reference genome.
#Therefore, we are uncertain about the quality of this available data.
#The mutations are from COSMIC and were associated by overlap with the enhancer sequences.

# 3. Coordinates and hg38
source("./scripts/bed2coord.R")
coord = bed2coord(CancerEnD[c(17,5:6)])
enh19 = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", coord)), CancerEnD[c(17,5:6)], coord, orig_assembly = "hg19")
enh19 = enh19[!duplicated(enh19),]

source("./scripts/hgtohg_onestep.R")
enh19to38 = hgtohg_onestep(enh19[2:5], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
tables_CancerEnD = hgtohg2tables("hg19", "hg38", enh19, enh19to38)
enh_CancerEnD = tables_CancerEnD$enhancers

# 4. Complete table
CancerEnD = cbind(CancerEnD[c(17,5,6)], orig_assembly = "hg19", CancerEnD[c(3,16)], enh_PMID = 24670763, 
                  biosample_name = NA, enh_method = "CAGE", type = "Transcribed enhancer", CancerEnD[c(15,2,1)], 
                  enh2gene_PMID = 32360910, enh2gene_method = "distance 0.1MB", hgnc_symbol_TFs = NA, TFs2enh_PMID = NA, 
                  TFs2enh_method = NA, CancerEnD[4], disease_PMID = 32360910, disease_method = NA, refseq_ID = NA, 
                  mutation_PMID = 32360910, mutation_method = NA, CancerEnD[7:12]) 
#disease_method -> NA and the mutations are in COSMIC format
CancerEnD = cbind(enh_ID = paste0("hg19_", gsub("[:-]", "_", coord)), CancerEnD)
CancerEnD = merge(enh_CancerEnD, CancerEnD[c(1,6:31)], by = "enh_ID")
CancerEnD$source = "CancerEnD"
indexes = which(CancerEnD$hasMutation %in% c("Confirmed somatic variant", "Variant of unknown origin"))
CancerEnD$mutation_method[indexes] = "overlap"

set = CancerEnD
set$enh2gene_method = "eQTL"
CancerEnD = rbind(CancerEnD, set)
rm(set, crossref, coord, indexes)

#5. Filtered table
filt_CancerEnD = CancerEnD[, 1:31]
filt_CancerEnD = filt_CancerEnD[,-19]
indexes = which(is.na(filt_CancerEnD$current_chr))
filt_CancerEnD = filt_CancerEnD[-indexes,]
filt_CancerEnD$mutation_PMID = filt_CancerEnD$mutation_method = NA
filt_CancerEnD[is.na(filt_CancerEnD)] = "-"

#Diseases
library(readxl)
diseases = read_excel("diseases.xlsx", col_names = T)
colnames(diseases)[1] = "disease"
library(tidyr)
diseases = pivot_longer(diseases[1:4], cols = c("OMIM", "DiseaseOntology", "MeSH"), names_to = "ontology", values_to = "ID")
indexes = which(diseases$ID == "-")
diseases = diseases[-indexes, -2]
filt_CancerEnD = merge(filt_CancerEnD, diseases, by = "disease", all.x = T)
filt_CancerEnD = filt_CancerEnD[c(2:25,31,26:30)]
colnames(filt_CancerEnD)[25] = "disease"
library(dplyr)
filt_CancerEnD = arrange(filt_CancerEnD, enh_ID)
rm(indexes)

#Diseases --> biosamples
colnames(diseases) = c("disease_name", "disease")
filt_CancerEnD = merge(filt_CancerEnD, diseases, by = "disease", all.x = T)
filt_CancerEnD = filt_CancerEnD[,-1]
colnames(filt_CancerEnD)[30] = "disease"
filt_CancerEnD = filt_CancerEnD[,c(1:24,30,25:29)]
filt_CancerEnD$biosample_name = filt_CancerEnD$disease
filt_CancerEnD$disease = "-"
filt_CancerEnD = filt_CancerEnD[!duplicated(filt_CancerEnD),]


#6. Save files
dir.create("./CancerEnD_results/")
save.image("./CancerEnD_results/CancerEnD_data.RData")
write.table(filt_CancerEnD, "./CancerEnD_results/CancerEnD.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
