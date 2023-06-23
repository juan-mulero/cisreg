#SCREEN V3 

library(data.table)
pELS = fread("https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.pELS.bed")
pELS$V6 = "pELS"
pELS = pELS[!duplicated(pELS),]

dELS = fread("https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.dELS.bed")
dELS$V6 = "dELS"
dELS = dELS[!duplicated(dELS),]

SCREEN = rbind(pELS, dELS)
source("./scripts/bed2coord.R")
coord = bed2coord(SCREEN[,1:3])
enh_ID = gsub("[:-]", "_", coord)
enh_ID = paste0("hg38_", enh_ID)
rm(coord)

crossref = paste0("https://screen.encodeproject.org/search/?q=", SCREEN$V5,"&assembly=GRCh38")
SCREEN = data.frame(enh_ID = enh_ID, orig_chr = SCREEN$V1, orig_start = SCREEN$V2, orig_end = SCREEN$V3, orig_assembly = "hg38",
                    current_chr = SCREEN$V1, current_start = SCREEN$V2, current_end = SCREEN$V3, current_assembly = "hg38",
                    minimum_ratio = "-", score = "-", original_ID = SCREEN$V5, crossref = crossref, enh_PMID = 32728249,
                    biosample_name = "-", enh_method = "DNase-seq,ChIP-seq", type = SCREEN$V6, source = "SCREEN",
                    hgnc_symbol_target_genes = "-", enh2gene_PMID = "-", enh2gene_method = "-", hgnc_symbol_TFs = "-",
                    TFs2enh_PMID = "-", TFs2enh_method = "-", disease = "-", disease_PMID = "-", disease_method = "-",
                    refsnp_ID = "-", mutation_PMID = "-", mutation_method = "-")

set = SCREEN
SCREEN$enh_method = "DNase-seq"
set$enh_method = "ChIP-seq"
SCREEN = rbind(SCREEN, set)
SCREEN = SCREEN[!duplicated(SCREEN),]
SCREEN = SCREEN[order(SCREEN$orig_chr, SCREEN$orig_start, SCREEN$orig_end),]
rm(crossref, enh_ID, set)

dir.create("SCREEN_results")
save.image("./SCREEN_results/SCREEN_data.RData")
write.table(SCREEN, "./SCREEN_results/SCREEN.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
