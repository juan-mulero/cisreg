#1. Reading and preprocess
library(rvest)
library(xml2)

#Download the data table
url = "https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?search.form=no;show=1;page_size=20000;form=search;action=search;search.result=yes;page=1"
read = read_html(url)
orig_VistaEnhancer = html_nodes(read, "table")
orig_VistaEnhancer = data.frame(html_table(orig_VistaEnhancer[6]))

VistaEnhancer = orig_VistaEnhancer
VistaEnhancer = VistaEnhancer[-1, 1:3]
colnames(VistaEnhancer) = c("original_ID", "coord", "external_gene_name")
indexes = which(grepl("^hs" ,VistaEnhancer$original_ID)) #To filter the mouse sequences
VistaEnhancer = VistaEnhancer[indexes,]
VistaEnhancer$external_gene_name = gsub("\\(intragenic\\)", "", VistaEnhancer$external_gene_name)

#The genes are separated by the "-" character. However, some genes incorporate this character
for (i in 1:length(VistaEnhancer$external_gene_name)){
  count = length(gregexpr("-", VistaEnhancer$external_gene_name[i])[[1]]) #Number of "-"
  if (count == 3){ #If there are 4 elements, we divide in the middle
    index = gregexpr("-", VistaEnhancer$external_gene_name[i])[[1]][2] #Position of the "-" to replace
    substr(VistaEnhancer$external_gene_name[i], index, index) = ","
  } else if(count == 1) { #If there are 2 elements, we divide in the middle
    VistaEnhancer$external_gene_name[i] = gsub("-", ",", VistaEnhancer$external_gene_name[i])
  }
  else { #If there are 3 elements, we need to preprocess
    indexes = gregexpr("-", VistaEnhancer$external_gene_name[i])[[1]] #Position of the "-"
    for (j in 1:length(indexes)){
      #We consider those situations in which genes have a dash. If this dash is not part of the gene name, it is replaced
      if (!(grepl("^(AS|IT|OT)", substr(VistaEnhancer$external_gene_name[i], indexes[j]+1, indexes[j]+2)) | #We check whether the word after the dash starts with AS (antisense), IT (intronic transcript) or OT (overlapping transcript)
            grepl("[0-9]", substr(VistaEnhancer$external_gene_name[i], indexes[j]+1, indexes[j]+1)) | #We check if the next character is a digit
            grepl("INS", substr(VistaEnhancer$external_gene_name[i], indexes[j]-3, indexes[j]-1)))){ #We check whether it is preceded by the characters INS (insulin)
        substr(VistaEnhancer$external_gene_name[i], indexes[j], indexes[j]) = ","
      }
    }
  }
}

#Coordinates
VistaEnhancer$coord = gsub(",", "", VistaEnhancer$coord)
source("./scripts/coord2bed.R")
coord = coord2bed(VistaEnhancer$coord)
VistaEnhancer = merge(coord, VistaEnhancer, by = "coord", all = T)
VistaEnhancer = VistaEnhancer[!duplicated(VistaEnhancer),]

#Crossref
ID = unlist(lapply(VistaEnhancer$original_ID, function(original_ID){unlist(strsplit(original_ID, "hs"))[2]}))
crossref = unlist(lapply(ID, function(id){
  paste0("https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?form=presentation&show=1&experiment_id=", id, "&organism_id=1")}))
VistaEnhancer = cbind(VistaEnhancer, crossref)
rm(coord, count, i, j, indexes, url, read, ID, crossref, index) 
  

#2. Enhancers:
enh19 = cbind(enh_ID = paste0("hg19_",gsub("[:-]", "_", VistaEnhancer$coord)), orig_chr = VistaEnhancer$chr, 
                  orig_start = VistaEnhancer$start, orig_end = VistaEnhancer$end, coord = VistaEnhancer$coord, 
                  orig_assembly = "hg19")
enh19 = data.frame(enh19[!duplicated(enh19),])

enh_metadata = cbind(enh_ID = paste0("hg19_",gsub("[:-]", "_", VistaEnhancer$coord)), VistaEnhancer[c(1,5,7)], 
                     enh_PMID = 17130149, enh_method = "sequence conservation,reporter gene assay", type = "enhancer", 
                     source = "VISTA")

#Biosamples:
orig_biosamples = read.delim("https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?form=search;page_size=20000;search.form=no;search.result=yes;action=search;page=1;show=1;search.sequence=1", 
                        header = F, sep = "\f")
biosamples = orig_biosamples
indexes = grep(">Human", biosamples$V1) #We only need the headers
biosamples = biosamples$V1[indexes]
biosamples = gsub("<pre>", "", biosamples)
coord = biosample_name = c()
count = 0
for (i in 1:length(biosamples)){ #We extract the positive biosamples from the header
  if (grepl("positive", biosamples[i])){
    count = count + 1
    split = unlist(strsplit(biosamples[i], "\\|"))
    coord[count] = trimws(split[2])
    split = split[5:length(split)]
    split = trimws(split)
    for (j in 1:length(split)){
      if (grepl("\\[", split[j])){
        split[j] = unlist(strsplit(split[j], "\\["))[1]
      }
    }
    biosample_name[count] = paste(split, collapse = ";") #We use ";" because the value trigeminal V (ganglion, cranial) have one comma
  }
}
data = data.frame(coord, biosample_name)
enh_metadata = merge(enh_metadata, data, by = "coord", all = T)
enh_metadata = enh_metadata[c(2:5,9,6:8)]
enh_metadata[is.na(enh_metadata)] = "-"
rm(biosample_name, biosamples, coord, count, i, j, indexes, split)


#Remap hg38
source("./scripts/hgtohg_onestep.R")
hg19tohg38 = hgtohg_onestep(enh19[2:5], 19, 38, 0.05)

source("./scripts/hgtohg2tables.R")
hg19tohg38_tables = hgtohg2tables("hg19", "hg38", enh19, hg19tohg38)


#Enhancer metadata split
enh_metadata_split = enh_metadata
source("./scripts/SplitIntoAtomicValues.R")
enh_metadata_split = SplitIntoAtomicValues(enh_metadata_split, 5, ";") #Biosamples
enh_metadata_split = SplitIntoAtomicValues(enh_metadata_split, 6, ",") #Methods
enh_metadata_split = enh_metadata_split[!duplicated(enh_metadata_split),]


#3. Target genes
VistaEnhancer = merge(VistaEnhancer, enh19[c(1,5)], by = "coord", all = T)
VistaEnhancer = cbind(VistaEnhancer, enh_PMID = 17130149, source = "VISTA")

enh2target_genes = cbind(VistaEnhancer[c(8,6)], enh2gene_method = "bracketing genes", VistaEnhancer[9:10])
enh2target_genes = SplitIntoAtomicValues(enh2target_genes, 2, ",")
enh2target_genes = enh2target_genes[!duplicated(enh2target_genes),]

upd_target_genes = unique(enh2target_genes$external_gene_name)
source("./scripts/UpdateGeneSymbol.R")
upd_target_genes = UpdateGeneSymbol(upd_target_genes)
enh2target_genes = merge(enh2target_genes, upd_target_genes, by = "external_gene_name", all = T)
enh2target_genes = merge(enh2target_genes, enh_metadata[c(1,5)], by = "enh_ID", all.x = T)
enh2target_genes = enh2target_genes[!duplicated(enh2target_genes),]
enh2target_genes = enh2target_genes[c(1:2,6,4,3,7,5)]
enh2target_genes[is.na(enh2target_genes)] = "-"
colnames(enh2target_genes)[4] = "enh2gene_PMID"

enh2target_genes_split = SplitIntoAtomicValues(enh2target_genes, 6, ";") #Biosamples
enh2target_genes_split = enh2target_genes_split[!duplicated(enh2target_genes_split),]


#4. Methods
method = c(enh_metadata_split$enh_method, enh2target_genes_split$enh2gene_method)
method = unique(method)
method = data.frame(method) #For the mapping


#5. Biosamples
biosample_name = c(enh_metadata_split$biosample_name, enh2target_genes_split$biosample_name)
biosample_name = unique(biosample_name)
biosample_name = data.frame(biosample_name)
index = which(biosample_name$biosample_name == "-")
if(length(index) > 0){biosample_name = biosample_name[-index,]}
biosample_name = data.frame(biosample_name) #For the mapping
rm(index)


#6. Merge data in an unique file
merge_VistaEnhancer = merge(hg19tohg38_tables$enhancers, enh_metadata_split, by = "enh_ID", all = T)
merge_VistaEnhancer = merge(merge_VistaEnhancer, enh2target_genes_split[1:5], by = "enh_ID", all = T)
merge_VistaEnhancer = merge_VistaEnhancer[!duplicated(merge_VistaEnhancer),]
colnames(merge_VistaEnhancer)[19:20] = c("external_gene_name_target_genes", "hgnc_symbol_target_genes")
merge_VistaEnhancer = cbind(merge_VistaEnhancer, external_gene_name_TFs = "-", hgnc_symbol_TFs = "-", TFs2enh_PMID = "-",
                            TFs2enh_method = "-", disease = "-", OMIM = "-", DO = "-", mesh = "-", disease_PMID = "-",
                            disease_method = "-", refsnp_ID = "-", mutation_PMID = "-", mutation_method = "-")
merge_VistaEnhancer[is.na(merge_VistaEnhancer)] = "-"

filt_merge_VistaEnhancer = merge_VistaEnhancer[-c(19,23,28:30)] #external_gene_names and diseases
filt_merge_VistaEnhancer = filt_merge_VistaEnhancer[!duplicated(filt_merge_VistaEnhancer),]
indexes = which(filt_merge_VistaEnhancer$hgnc_symbol_target_genes == "-")
if (length(indexes) > 0){filt_merge_VistaEnhancer$enh2gene_PMID[indexes] = filt_merge_VistaEnhancer$enh2gene_method[indexes] = "-"}
rm(indexes)


#7. Save
dir.create("./VistaEnhancer_results")
save.image("./VistaEnhancer_results/Data_VistaEnhancer.RData")
write.table(filt_merge_VistaEnhancer, "./VistaEnhancer_results/VISTA.tsv", quote = F, sep = "\t", 
            col.names = T, row.names = F, fileEncoding = "utf8")
