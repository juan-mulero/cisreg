#Data table of RAEdb in download section --> http://www.computationalbiology.cn/RAEdb/download.html
#PMID Study Method  Species CellLine    GEO Accession
#23328393	Arnold CD,2015	STARR-seq	Human	HeLa cells      GSE40739		
#26367798	Rathert P,2015	STARR-seq	Human	K562 cells      GSE63782		
#29151363	Liu Y,2017	STARR-seq	Human	LNCaP cells     GSE82204		
#29256496	Felix Muerdter,2017	STARR-seq	Human	HeLaS3 cells    GSE100432		
#27701403	Ernst J,2016	MPRA	Human	HepG2 cells	GSE71279	
#27701403	Ernst J,2016	MPRA	Human	K562 cells	GSE71279	
#27524623	Groff AF,2016	MPRA	Mouse	C2C12 cells	GSE73472	
#28973438	Younger ST,2017	MPRA	Human	Fibroblasts	GSE83780	
#27831498	Inoue F,2017	MPRA	Human	HepG2 cells	GSE83894	

#URL = c("http://www.computationalbiology.cn/RAEdb/download/GSE40739_HeLa_bed.rar",
#        "http://www.computationalbiology.cn/RAEdb/download/GSE63782_K562_bed.rar",
#        "http://www.computationalbiology.cn/RAEdb/download/GSE82204_LNCaP_bed.rar",
#        "http://www.computationalbiology.cn/RAEdb/download/GSE100432_HeLaS3_bed.rar",
#        "http://www.computationalbiology.cn/RAEdb/download/GSE71279_bed.rar",
#        "http://www.computationalbiology.cn/RAEdb/download/GSE71279_bed.rar",
#        "http://www.computationalbiology.cn/RAEdb/download/GSE83780_Fibroblasts_bed.rar",
#        "http://www.computationalbiology.cn/RAEdb/download/GSE83894_HepG2_bed.rar")

#Reading of files
library(data.table)
enh_PMID = c(23328393, 26367798, 29151363, 29256496, 27701403, 27701403, 28973438, 27831498)
biosample_name = c("HeLa", "K562", "LNCaP", "HeLaS3", "HepG2", "K562", "Fibroblasts", "HepG2")
enh_method = c("STARR-seq", "STARR-seq", "STARR-seq", "STARR-seq", "MPRA", "MPRA", "MPRA", "MPRA")
GEO = c("GSE40739", "GSE63782", "GSE82204", "GSE100432", "GSE71279", "GSE71279", "GSE83780", "GSE83894")
df_metadata = data.frame(enh_PMID, biosample_name, enh_method, GEO)
list_dirs = list.dirs("./RAEdb_data")
index = which(list_dirs == "./RAEdb_data")
if (length(index) > 0){list_dirs = list_dirs[-index]}
df_enh = df_epr = data.frame()
for (i in 1:length(list_dirs)){
        dir = list_dirs[i]
        dir_split = unlist(strsplit(dir, "/|_"))
        GEO = dir_split[grep("GSE", dir_split)]
        list_files = list.files(dir)
        for (j in 1:length(list_files)){
                file = list_files[j]
                file_split = unlist(strsplit(file, "_|\\."))
                biosample_name = file_split[1]
                type = file_split[length(file_split)-1]
                data = data.frame(fread(paste0(dir,"/",file)))
                if (type == "enhancer"){df_enh = rbind(df_enh, cbind(data[1:4], biosample_name, type, GEO))}
                else if (type == "epromoter"){df_epr = rbind(df_epr, cbind(data[1:3], biosample_name, type, GEO))}
                else {print("Error, the file is not in the correct format")}
        }
}

#Unique table
df_enh = df_enh[!duplicated(df_enh),]
colnames(df_enh)[1:4] = c("chr", "start", "end", "activity")
df_epr = df_epr[!duplicated(df_epr),]
colnames(df_epr)[1:3] = c("chr", "start", "end")
df_epr = cbind(df_epr[1:3], activity = NA, df_epr[4:6])

cisreg = rbind(df_enh, df_epr)
cisreg = merge(cisreg, df_metadata, by = c("biosample_name", "GEO"), all.x = T)
cisreg = cisreg[c(3:5,8,1,9,7,2)]

#Enrich table
source("./scripts/bed2coord.R")
coord = bed2coord(cisreg[1:3])
RAEdb = cbind(enh_ID = paste0("hg38_", gsub("[:-]", "_", coord)), cisreg[1:3], orig_assembly = "hg38", cisreg[1:3], 
              current_assembly = "hg38", minimum_ratio = NA, score = NA, original_ID = NA, crossref = NA, cisreg[4:7],
              source = "RAEdb", hgnc_symbol_taget_genes = NA, enh2gene_PMID = NA, enh2gene_method = NA, hgnc_symbol_TFs = NA,
              TFs2enh_PMID = NA, TFs2enh_method = NA, disease = NA, disease_PMID = NA, disease_method = NA, refseq_ID = NA,
              mutation_PMID = NA, mutation_method = NA)
colnames(RAEdb)[c(2:4,6:8)] = c("orig_chr", "orig_start", "orig_end", "current_chr", "current_start", "current_end")
RAEdb[is.na(RAEdb)] = "-"

#Save and write
dir.create("./RAEdb_results")
rm(biosample_name, coord, dir, dir_split, enh_method, enh_PMID, file, file_split, GEO, i, index, j, list_dirs, list_files, type, data)
save.image("./RAEdb_results/data_RAEdb.RData")
write.table(RAEdb, "./RAEdb_results/RAEdb.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
