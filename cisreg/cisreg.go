package main

import(
        "fmt"
        "bufio"
        "os"
	"io/ioutil"
        "strings"
	//log "github.com/sirupsen/logrus"
        "github.com/vlmir/bgw3/src/semweb"
        "github.com/vlmir/bgw3/src/util"
	"github.com/vlmir/bgw3/src/bgw"
	//"encoding/json"
	//"reflect"
	"strconv"
)

//Set of human chromosomes
var chromosomes = util.SliceSet{
        "chr1": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000001.11>", "chr-1", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr2": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000002.12>", "chr-2", "<http://purl.obolibrary.org/obo/SO_0000340>"},
	"chr3": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000003.12>", "chr-3", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr4": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000004.12>", "chr-4", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr5": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000005.10>", "chr-5", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr6": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000006.12>", "chr-6", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr7": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000007.14>", "chr-7", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr8": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000008.11>", "chr-8", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr9": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000009.12>", "chr-9", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr10":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000010.11>", "chr-10", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr11":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000011.10>", "chr-11", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr12":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000012.12>", "chr-12", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr13":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000013.11>", "chr-13", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr14":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000014.9>", "chr-14", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr15":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000015.10>", "chr-15", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr16":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000016.10>", "chr-16", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr17":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000017.11>", "chr-17", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr18":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000018.10>", "chr-18", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr19":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000019.10>", "chr-19", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr20":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000020.11>", "chr-20", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr21":        {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000021.9>", "chr-21", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chr22":	{"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000022.11>", "chr-22", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chrX": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000023.11>", "chr-X", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chrY": {"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000024.10>", "chr-Y", "<http://purl.obolibrary.org/obo/SO_0000340>"},
}

//Map of diseases to URI
var diseases2uri = map[string]string{
	"DOID": rdf.Nss["doid"],
        "OMIM": rdf.Nss["omim"],
        "mesh": rdf.Nss["mesh"],
}

//Map of sources to URI
var sources2uri = util.SliceSet{
        "ENdb": {"<http://www.licpathway.net/ENdb/>", "ENdb", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
        "EnDisease":    {"<http://health.tsinghua.edu.cn/jianglab/endisease/>", "EnDisease", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
        "FANTOM5":      {"<https://fantom.gsc.riken.jp/5/>", "FANTOM5", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
        "VISTA":        {"<https://enhancer.lbl.gov/>", "VISTA", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"dbSUPER":      {"<https://asntech.org/dbsuper/>", "dbSUPER", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"Ensembl":	{"<https://www.ensembl.org/Homo_sapiens/Info/Index>", "Ensembl", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"RefSeq":	{"<https://www.ncbi.nlm.nih.gov/refseq/>", "RefSeq", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"RAEdb":	{"<http://www.computationalbiology.cn/RAEdb/index.php>", "RAEdb", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"EnhancerDB":	{"<http://lcbb.swjtu.edu.cn/EnhancerDB/>", "EnhancerDB", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"ChromHMM":	{"<https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeBroadHmm&db=hg19>", "ChromHMM", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"TADKB":	{"<http://dna.cs.miami.edu/TADKB/>", "TADKB", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"3DGB":	{"<http://3dgenome.fsm.northwestern.edu/index.html>", "3DGB", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
}

//Map of assembly versions
var assembly2uri = util.SliceSet{
        "hg38": {"<https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39>", "GRCh38.p13", "<http://purl.obolibrary.org/obo/SO_0001248>"},
	//"hg19": {"<https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25>", "GRCh37.p13", "<http://purl.obolibrary.org/obo/SO_0001248>"},
}

//Datatype properties
var Dpys = util.SliceSet{
	"seq2start":    {"<http://purl.obolibrary.org/obo/OGI_1000004>", "start point of interval"},
	"seq2end":      {"<http://purl.obolibrary.org/obo/OGI_1000003>", "end point of interval"},
	"seq2strand":  {"<http://purl.obolibrary.org/obo/GENO_0000906>", "on strand"},
}

//BGW properties to Biolink properties
var propbgw2propbiolink = util.SliceSet{
	"sth2lbl":	{"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/linkml/alias>"},
	"sub2cls":	{"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/linkml/is_a>"},
	"sub2ppy":	{"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/linkml/is_a>"},
	"sth2mtd":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/has_evidence>"},
	"sth2ori":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/download_url>"},
	"seq2version":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/genome_build>"},
	"sth2evd":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/publication_id>"},
	"sth2src": {"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/supporting_data_source>"},
	"rgr2trg":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/genotype_to_gene_association_object>"},
	"tlp2tlp":	{"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/biolink/vocab/directly_physically_interacts_with>"},
	"gn2phn":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/genotype_to_disease_association_object>"},
	"gn2txn":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/in_taxon>"},
	"seq2chr":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/genomic_sequence_localization_object>"},
	"seq2cell_line":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/supporting_study_context>"},
	"seq2cell":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/supporting_study_context>"},
	"seq2anato":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/supporting_study_context>"},
	"seq2start":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/start_coordinate>"},
	"seq2end":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/end_coordinate>"},
	"seq2strand":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/strand>"},
}

//Func to check if one number is present in one slice
func ExistInSlice(set []int, search int) bool{
        if len(set) > 0 {
		for _,value := range set{
                	if value == search{
                        	return true
                	}
        	}
		return false
	}
        return false
}

//Enhancer table to map3D
	//pathfile : path of the file to convert (input)
	//util.Set3D : map3D with de data obtained (output)
func Enh2map3D(pathfile string, p_old2new_IDs, p_new2old_IDs *util.Set3D) (enhmap, old2new_IDs, new2old_IDs util.Set3D, err error){
        file, err := os.Open(pathfile)
        if err != nil{
                panic(err)
	}

	old2new_IDs = *p_old2new_IDs
	new2old_IDs = *p_new2old_IDs
	enhmap = make(util.Set3D)
        name_fields := []string{        //column names to names that we will use like secondary keys in the Set3D
                        "old_ID",
                        "orig_chr",
                        "orig_start",
                        "orig_end",
                        "orig_assembly",
                        "chr",
                        "start",
                        "end",
                        "assembly",
                        "min_ratio",
                        "score",
                        "orig_id",
                        "crossref",
                        "pubmed",
                        "biosample_name",
                        "method",
                        "type",
                        "source",
                        "target_gene",
                        "pubmed",
                        "method",
                        "tfactor",
                        "pubmed",
                        "method",
                        "disease",
                        "pubmed",
                        "method",
                        "refsnp",
                        "pubmed",
                        "method",
                        }

        skip := []int{0,1,2,3,4,9,10,11,16,26,27,28,29}         //Columns that we will not implement at this moment

	scanner := bufio.NewScanner(file)
        i := 0

        for scanner.Scan() {
                i++
                line := scanner.Text()
                if line == "" {
                        panic(fmt.Sprintf("Empty line in the row number %v", i))
                }

                if i != 1 {     //Because the files have header with the names of the different columns
                        fields := strings.Split(line, "\t")
			chr := strings.Split(fields[5], "chr")[1]
			assembly := strings.Split(fields[8], "hg")[1]
			fields_enhID := []string{assembly, chr, fields[6], fields[7]}    //Fields that contain the values to construct the ID
                        enh_ID := strings.Join(fields_enhID, ".")   //Enhancer ID as current coordinates

	                var primarykey string
        	        if _,ok := old2new_IDs[enh_ID]; ok {
                	        primarykey = old2new_IDs[enh_ID]["new_ID"].Keys()[0]
	                } else {
        	                n_ids := len(new2old_IDs.Keys()) + 1
                	        n_zeros := 11 - len(strconv.Itoa(n_ids))
                        	primarykey = strings.Join([]string{"CRMHS", strings.Repeat("0", n_zeros), strconv.Itoa(n_ids)}, "")
                        	old2new_IDs.Add(enh_ID, "new_ID", primarykey)
                        	new2old_IDs.Add(primarykey, "old_ID", enh_ID)
                	}

			label := strings.Join(fields[6:8], "-")
			label = strings.Join([]string{`"`, fields[5], ":", label, `"`}, "")
			enhmap.Add(primarykey, "label", label)

                        for j,v := range name_fields{
                                if !ExistInSlice(skip, j) {
                                        secondarykey := v
                                        tertiarykey := fields[j]
                                        enhmap.Add(primarykey, secondarykey, tertiarykey)
                                }
                        }
                }
        }
        return enhmap, old2new_IDs, new2old_IDs, nil
}

//TAD table to map3D
        //pathfile : path of the file to convert (input)
        //util.Set3D : map3D with de data obtained (output)
func TAD2map3D(pathfile string, p_old2new_IDs, p_new2old_IDs *util.Set3D) (tadmap, old2new_IDs, new2old_IDs util.Set3D, err error){
        file, err := os.Open(pathfile)
        if err != nil{
                panic(err)
        }

        old2new_IDs = *p_old2new_IDs
        new2old_IDs = *p_new2old_IDs
        tadmap = make(util.Set3D)
        name_fields := []string{        //column names to names that we will use like secondary keys in the Set3D
                        "old_ID",
                        "orig_chr",
                        "orig_start",
                        "orig_end",
                        "orig_assembly",
                        "chr",
                        "start",
                        "end",
                        "assembly",
                        "min_ratio",
                        "score",
                        "orig_id",
                        "crossref",
                        "pubmed",
                        "biosample_name",
                        "method",
			"resolution",
                        "type",
                        "source",
			}
	skip := []int{0,1,2,3,4,9,10,11,16,17}         //Columns that we will not implement at this moment

        scanner := bufio.NewScanner(file)
        i := 0

        for scanner.Scan() {
                i++
                line := scanner.Text()
                if line == "" {
                        panic(fmt.Sprintf("Empty line in the row number %v", i))
                }

                if i != 1 {     //Because the files have header with the names of the different columns
                        fields := strings.Split(line, "\t")
                        chr := strings.Split(fields[5], "chr")[1]
                        assembly := strings.Split(fields[8], "hg")[1]
                        fields_tadid := []string{assembly, chr, fields[6], fields[7]}    //Fields that contain the values to construct the ID
                        tad_id := strings.Join(fields_tadid, ".")   //TAD ID as current coord

                        var primarykey string
                        if _,ok := old2new_IDs[tad_id]; ok {
                                primarykey = old2new_IDs[tad_id]["new_ID"].Keys()[0]
                        } else {
                                n_ids := len(new2old_IDs.Keys()) + 1
                                n_zeros := 11 - len(strconv.Itoa(n_ids))
                                primarykey = strings.Join([]string{"TADHS", strings.Repeat("0", n_zeros), strconv.Itoa(n_ids)}, "")
                                old2new_IDs.Add(tad_id, "new_ID", primarykey)
                                new2old_IDs.Add(primarykey, "old_ID", tad_id)
                        }

                        label := strings.Join(fields[6:8], "-")
                        label = strings.Join([]string{`"`, fields[5], ":", label, `"`}, "")
                        tadmap.Add(primarykey, "label", label)

			for j,v := range name_fields{
                                if !ExistInSlice(skip, j) {
                                        secondarykey := v
                                        tertiarykey := fields[j]
                                        tadmap.Add(primarykey, secondarykey, tertiarykey)
                                }
                        }
                }
        }
        return tadmap, old2new_IDs, new2old_IDs, nil
}

//Chromosomes table to map
        //pathfile : path of the file to convert (input)
        //util.Set4D : map with de data obtained (output)
func Chr2map(pathfile string) (util.Set4D, error){
        file, err := os.Open(pathfile)
        if err != nil{
                panic(err)
        }

        map4D := make(util.Set4D)
        name_fields := []string{        //column names
                        "txid",
                        "chr",
                        "chr_label",
                        "chr_urirefseq",
			"assembly_label",
			"assembly_urirefseq",
                        }
        skip := []int{0,1}         //Columns to skip

        scanner := bufio.NewScanner(file)
        i := 0
	for scanner.Scan() {
		i++
                line := scanner.Text()
                if line == "" {
                        panic(fmt.Sprintf("Empty line in the row number %v", i))
                }

                fields := strings.Split(line, "\t")
                primarykey := fields[0]
		secondarykey := fields[1]

                for j,v := range name_fields{
                        if !ExistInSlice(skip, j){
                                tertiarykey := v
                                quaternarykey := fields[j]
                                map4D.Add(primarykey, secondarykey, tertiarykey, quaternarykey)
                        }
                }
        }
        return map4D, nil
}

//Func to convert one file with enrichment information to map 3D
	//path : path of the file to convert (input)
	//name_fields : name of the columns of the table that will be use like secondary keys (input)
	//skip : slice with the number of the columns that we do not want to use (start in 0) (input)
	//util.Set3d : map3D with the data obtained (output)
func EnrichFile2map3D(path string, name_fields []string, skip []int) (util.Set3D, error){
        file, err := os.Open(path)
        if err != nil{
                panic(err)
        }
        map3D := make(util.Set3D)
        scanner := bufio.NewScanner(file)
        i := 0
        for scanner.Scan() {
                i++
                line := scanner.Text()
                if line == "" {
                        panic(fmt.Sprintf("Empty line in the row number %v", i))
                }

                fields := strings.Split(line, "\t")
                primarykey := fields[0]
                for j,v := range name_fields{
                        if j!= 0 && !ExistInSlice(skip, j){
                                secondarykey := v
                                tertiarykey := fields[j]
                                map3D.Add(primarykey, secondarykey, tertiarykey)
                        }
                }
        }
        return map3D, nil
}

//Func to export the table of the enhancers to an RDF file (use the set3D as intermediary)
	//path_source : path of the file, the table with the enhancers (input)
	//path_rdf : path where we want to save the RDF file (input)
	//name_source : name of the database used as source
	//map_triplets : map3D with the triplets generated for the RDF file
	//old2new_IDs: map3D with the old IDs and the new IDs
	//new2old_IDs: map3D with the new and the old IDs
	//num_ln : number of triplets of the RDF file
func ExportEnh2rdf(path_source, path_rdf, source_name string, p_genes_xmap *bgw.Xmap, p_old2new_IDs, p_new2old_IDs *util.Set3D)	(map_triplets, old2new_IDs, new2old_IDs util.Set3D, num_ln int, err error){
	//File creation
        rdfFile2export, err := os.Create(path_rdf)
        if err != nil {
                panic(err)
        }
        defer rdfFile2export.Close()
        var output_file strings.Builder
        keys4cisreg := make(util.SliceSet)

        //Annotation properties
        keys4cisreg["Apys"] = []string{
                "sth2dfn",
                "sth2lbl",
                "sth2syn",
        }

        //Object properties
        keys4cisreg["Opys"] = []string{
                "ins2cls",
                "sub2cls",
                "sth2ori",
                "sth2evd",
                "sth2mtd",
                "sth2src",
                "reg2targ",
                "tlp2tlp",
                "gn2phn",
		"gn2txn",
		"seq2chr",
                "seq2version",
		"seq2cell_line",
		"seq2cell",
		"seq2anato",
		"biolinkcat",
		"is_a",
		"sth2exm",
        }

        //Parental classes
        keys4cisreg["Prns"] = []string{
                "crm",		//cis-regulatory module
		"tfactor",	//Transcription factor
		"assembly",	//genome assembly
		"chr",		//chromosomes
		"db",		//database
        }

        //Header of the RDF file with the definition of the properties and the parental classes
        num_ln = 0
        header, num_ln_header := rdf.Capita(keys4cisreg)
        output_file.WriteString(header)
        num_ln += num_ln_header

	//URIs for cisreg
        uris_cisreg := rdf.FmtURIs(keys4cisreg)

	//BGW properties to Biolink properties
	proptypes := keys4cisreg.Keys()
	for _,proptype := range proptypes {
		if proptype == "Apys" || proptype == "Opys" {
			properties := keys4cisreg[proptype]
			for _,property := range properties {
				if _, ok := propbgw2propbiolink[property]; ok {
					triplet := rdf.FormT(uris_cisreg[property], propbgw2propbiolink[property][0], propbgw2propbiolink[property][1])
					output_file.WriteString(triplet)
					num_ln++
				}
			}
		}
	}

	//Datatype properties (the function rdf.Capita does not include this type of properties)
	datatype_prop := []string{
		"seq2start",
		"seq2end",
	}
	for _,value := range datatype_prop {
		triplet := rdf.FormT(Dpys[value][0], uris_cisreg["ins2cls"], rdf.CompU(rdf.Nss["owl"], "DatatypeProperty"))
		output_file.WriteString(triplet)

		triplet = rdf.FormT(Dpys[value][0], uris_cisreg["sth2lbl"], rdf.FormL(Dpys[value][1]))
		output_file.WriteString(triplet)

		triplet = rdf.FormT(Dpys[value][0], propbgw2propbiolink[value][0], propbgw2propbiolink[value][1])
                output_file.WriteString(triplet)
                num_ln = num_ln + 3
	}

        //AltLabel CRM and biolink category
        triplet := rdf.FormT(uris_cisreg["crm"], uris_cisreg["sth2syn"], rdf.FormL("cis-regulatory module"))
        output_file.WriteString(triplet)

	triplet = rdf.FormT(uris_cisreg["crm"], uris_cisreg["is_a"], rdf.CompU(rdf.Nss["biolink"], "NucleicAcidEntity"))
        output_file.WriteString(triplet)
        num_ln = num_ln + 2

	//AltLabel TF
        triplet = rdf.FormT(uris_cisreg["tfactor"], uris_cisreg["sth2syn"], rdf.FormL("TF"))
        output_file.WriteString(triplet)
	num_ln++

	//Chromosomes
	keys_dic := chromosomes.Keys()
	for _,value := range keys_dic {
        	triplet := rdf.FormT(chromosomes[value][0], uris_cisreg["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
                output_file.WriteString(triplet)

                triplet = rdf.FormT(chromosomes[value][0], uris_cisreg["sub2cls"], chromosomes[value][2])
                output_file.WriteString(triplet)

                triplet = rdf.FormT(chromosomes[value][0], uris_cisreg["sth2lbl"], rdf.FormL(chromosomes[value][1]))
                output_file.WriteString(triplet)

		triplet = rdf.FormT(chromosomes[value][0], uris_cisreg["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "GenomicSequenceLocalization"))
                output_file.WriteString(triplet)
                num_ln = num_ln + 4
	}

	//Assemblies
	keys_dic = assembly2uri.Keys()
        for _,value := range keys_dic {
                triplet := rdf.FormT(assembly2uri[value][0], uris_cisreg["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
                output_file.WriteString(triplet)

                triplet = rdf.FormT(assembly2uri[value][0], uris_cisreg["sub2cls"], assembly2uri[value][2])
                output_file.WriteString(triplet)

                triplet = rdf.FormT(assembly2uri[value][0], uris_cisreg["sth2lbl"], rdf.FormL(assembly2uri[value][1]))
                output_file.WriteString(triplet)

		triplet = rdf.FormT(assembly2uri[value][0], uris_cisreg["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "Genome"))
		output_file.WriteString(triplet)
                num_ln = num_ln + 4
        }

	//Database
	triplet = rdf.FormT(sources2uri[source_name][0], uris_cisreg["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
        output_file.WriteString(triplet)

	triplet = rdf.FormT(sources2uri[source_name][0], uris_cisreg["sub2cls"], sources2uri[source_name][2])
        output_file.WriteString(triplet)

	triplet = rdf.FormT(sources2uri[source_name][0], uris_cisreg["sth2lbl"], rdf.FormL(sources2uri[source_name][1]))
	output_file.WriteString(triplet)

        triplet = rdf.FormT(sources2uri[source_name][0], uris_cisreg["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "InformationResource"))
	output_file.WriteString(triplet)
        num_ln = num_ln + 4

        //source to RDF
        resource, old2new_IDs, new2old_IDs, err := Enh2map3D(path_source, p_old2new_IDs, p_new2old_IDs)
        if err != nil {
                panic(err)
        }

        //Map secondary keys of source to URI properties
        seckey2uri := map[string]string {
                "chr":		uris_cisreg["seq2chr"],
                "start":	Dpys["seq2start"][0],
                "end":		Dpys["seq2end"][0],
                "assembly":	uris_cisreg["seq2version"],
		"label":	uris_cisreg["sth2lbl"],
                //"orig_id":	uris_cisreg["sth2syn"],
                "crossref":	uris_cisreg["sth2ori"],
                "pubmed":	uris_cisreg["sth2evd"],
                "method":	uris_cisreg["sth2mtd"],
                "source":	uris_cisreg["sth2src"],
                "target_gene":	uris_cisreg["reg2targ"],
                "tfactor":	uris_cisreg["tlp2tlp"],
                "disease":	uris_cisreg["gn2phn"],
        }

        //Map 3D of methods
	name_fields := []string{
		"method_label",
		"efo",
		"obi",
		"eco",
		"bao",
		"ncit",
		"mi",
	}
	var skip []int
        methods2ids,err := EnrichFile2map3D("/home/jmulero/data/enrich/methods.tsv", name_fields, skip)
        if err != nil {
                panic(err)
        }

        //Map of method to URI properties
        methkey2uri := map[string]string{
                "efo":  rdf.Nss["efo"],
                "obi":  rdf.Nss["obo"],
                "eco":  rdf.Nss["obo"],
                "bao":  rdf.Nss["bao"],
                "ncit": rdf.Nss["obo"],
                "mi":   rdf.Nss["obo"],
        }

        //Map for gene enrich
	name_fields = []string{
		"symbol",
		"ncbigene",
                "uniprot",
	}
        symbol2enrich, err := EnrichFile2map3D("/home/jmulero/data/enrich/symbol2enrich.tsv", name_fields, skip)
        if err != nil {
                panic(err)
        }

	//Map for biosamples
	name_fields = []string{
		"label",
		"CLO",
		"CL",
		"UBERON",
		"BTO",
		"type",
	}
	biosamples, err := EnrichFile2map3D("/home/jmulero/data/enrich/biosamples.tsv", name_fields, skip)
	if err != nil {
                panic(err)
        }

        //Set3D of the source to triplets
	map_triplets = make(util.Set3D)
	str_cisreg9606 := "http://rdf.biogateway.eu/crm/9606/"
	genes_xmap := *p_genes_xmap
	uriref_cls := rdf.CompU(rdf.Nss["owl"], "Class")
	enh_ids := resource.Keys()
        for _,primarykey := range enh_ids {
		set2D := resource[primarykey]

                //The enh is a class
                enhid_uriref := rdf.CompU(str_cisreg9606, primarykey)
                triplet := rdf.FormT(enhid_uriref, uris_cisreg["ins2cls"], uriref_cls)
		map_triplets.Add(enhid_uriref, uris_cisreg["ins2cls"], uriref_cls)
                output_file.WriteString(triplet)
                num_ln++

                //The enh is a subclass of CRM
                triplet = rdf.FormT(enhid_uriref, uris_cisreg["sub2cls"], uris_cisreg["crm"])
		map_triplets.Add(enhid_uriref, uris_cisreg["sub2cls"], uris_cisreg["crm"])
                output_file.WriteString(triplet)
                num_ln++

		//Definition of the class
		var chr, start, end string
		for label,_ := range set2D["label"]{
			label = strings.Replace(label, "\"", "", -1)
			split := strings.Split(label, ":")
			chr = split[0]
			split = strings.Split(split[1], "-")
			start = split[0]
			end = split[1]
		}
		definition := fmt.Sprintf("Cis-regulatory module located in Homo sapiens %s between %s and %s", chr, start, end)
		triplet = rdf.FormT(enhid_uriref, uris_cisreg["sth2dfn"], rdf.FormL(definition))
		map_triplets.Add(enhid_uriref, uris_cisreg["sth2dfn"], rdf.FormL(definition))
                output_file.WriteString(triplet)
                num_ln++

		//prefLabel of classes
		triplet = rdf.FormT(enhid_uriref, uris_cisreg["sth2lbl"], rdf.FormL(strings.Join([]string{"crm/", primarykey}, "")))
                map_triplets.Add(enhid_uriref, uris_cisreg["sth2lbl"], rdf.FormL(strings.Join([]string{"crm/", primarykey}, "")))
                output_file.WriteString(triplet)
                num_ln++

		//Instances
                var instid_uriref, database string
                for source,_ := range set2D["source"]{
			database = strings.ToLower(source)
                        instid_uriref = fmt.Sprintf("<%s%s#%s>", str_cisreg9606, primarykey, database)
                }
                map_triplets.Add(instid_uriref, uris_cisreg["ins2cls"], enhid_uriref)
                inst2class := rdf.FormT(instid_uriref, uris_cisreg["ins2cls"], enhid_uriref)
                set_inst := fmt.Sprintf("%s", inst2class)
		num_ln++

		//prefLabel of instances
		triplet = rdf.FormT(instid_uriref, uris_cisreg["sth2lbl"], rdf.FormL(strings.Join([]string{"crm#", primarykey, "#", database}, "")))
                map_triplets.Add(instid_uriref, uris_cisreg["sth2lbl"], rdf.FormL(strings.Join([]string{"crm#", primarykey, "#", database}, "")))
		set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

		//The sequences belong to Homo sapiens
		humantaxon_uri := rdf.CompU(rdf.Nss["obo"], "NCBITaxon_9606")
		triplet = rdf.FormT(enhid_uriref, uris_cisreg["gn2txn"], humantaxon_uri)
	        map_triplets.Add(enhid_uriref, uris_cisreg["gn2txn"], humantaxon_uri)
                output_file.WriteString(triplet)

		triplet = rdf.FormT(instid_uriref, uris_cisreg["gn2txn"], humantaxon_uri)
		map_triplets.Add(instid_uriref, uris_cisreg["gn2txn"], humantaxon_uri)
                set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
		num_ln = num_ln + 2

		//NCBITaxon is a category of biolink:OrganismTaxon
		biolink_taxon := rdf.CompU(rdf.Nss["biolink"], "OrganismTaxon")
		triplet = rdf.FormT(humantaxon_uri, uris_cisreg["biolinkcat"], biolink_taxon)
		map_triplets.Add(humantaxon_uri, uris_cisreg["biolinkcat"], biolink_taxon)
		if map_triplets[humantaxon_uri][uris_cisreg["biolinkcat"]][biolink_taxon] == 1 {
			output_file.WriteString(triplet)
	                num_ln++
		}

                //Loop for secondary and tertiary keys
		secondarykeys := set2D.Keys()
                for _,secondarykey := range secondarykeys {
			set1D := set2D[secondarykey]
			tertiarykeys := set1D.Keys()
                        for _,tertiarykey := range tertiarykeys {
                                object, biolinkobject := "-", "-"
                                if  tertiarykey != "-" {
                                        if secondarykey == "pubmed" {
                                                object = rdf.CompU(rdf.Nss["pubmed"], tertiarykey)
						triplet := rdf.FormT(enhid_uriref, seckey2uri[secondarykey], object)
						triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(enhid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                output_file.WriteString(triplet)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln = num_ln + 2

                                                old_object := object
                                                object = rdf.CompU("http://www.ncbi.nlm.nih.gov/pubmed/", tertiarykey)
                                                triplet = rdf.FormT(old_object, uris_cisreg["sth2exm"], object)
                                                map_triplets.Add(old_object, uris_cisreg["sth2exm"], object)
                                                if map_triplets[old_object][uris_cisreg["sth2exm"]][object] == 1 {
                                                        output_file.WriteString(triplet)
                                                        num_ln++
                                                }

						biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Article")
						biolinktriplet := rdf.FormT(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                map_triplets.Add(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                if map_triplets[object][uris_cisreg["biolinkcat"]][biolinkobject] == 1 {
                                                	output_file.WriteString(biolinktriplet)
                                                        num_ln++
                                                }
						continue
                                        } else if secondarykey == "disease" {
						split := strings.Split(tertiarykey, ":")
						if uri, ok := diseases2uri[split[0]]; ok {
							object = rdf.CompU(uri, split[1])
							triplet := rdf.FormT(enhid_uriref, seckey2uri[secondarykey], object)
							triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                        map_triplets.Add(enhid_uriref, seckey2uri[secondarykey], object)
                                                        map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
							output_file.WriteString(triplet)
                                                        set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                        num_ln = num_ln + 2

							if split[0] == "mesh" {
								new_object := rdf.CompU("http://identifiers.org/mesh/", split[1])
								triplet = rdf.FormT(object, uris_cisreg["sth2exm"], new_object)
								map_triplets.Add(object, uris_cisreg["sth2exm"], new_object)
								if map_triplets[object][uris_cisreg["sth2exm"]][new_object] == 1 {
									output_file.WriteString(triplet)
                                                                	num_ln++
								}
							} else if split[0] == "OMIM" {
                                                                new_object := rdf.CompU("http://purl.obolibrary.org/obo/OMIM_", split[1])
                                                                triplet = rdf.FormT(object, uris_cisreg["sth2exm"], new_object)
                                                                map_triplets.Add(object, uris_cisreg["sth2exm"], new_object)
                                                                if map_triplets[object][uris_cisreg["sth2exm"]][new_object] == 1 {
                                                                        output_file.WriteString(triplet)
                                                                        num_ln++
                                                                }
								object = new_object
							}
                                                        biolinkobject := rdf.CompU(rdf.Nss["biolink"], "Disease")
                                                        biolinktriplet := rdf.FormT(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                        map_triplets.Add(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                        if map_triplets[object][uris_cisreg["biolinkcat"]][biolinkobject] == 1 {
                                                                output_file.WriteString(biolinktriplet)
                                                                num_ln++
                                                        }
						}
						continue
                                        } else if secondarykey == "method" {
                                                if _, ok := methods2ids[tertiarykey]; ok {
							ontologies := methods2ids[tertiarykey].Keys()
                                                        for _,ontology := range ontologies {
								methods_ids := methods2ids[tertiarykey][ontology].Keys()
                                                                for _,method_id := range methods_ids {
                                                                        if method_id != "-" {
										object = rdf.CompU(methkey2uri[ontology], method_id)
                                                                                triplet := rdf.FormT(enhid_uriref, seckey2uri[secondarykey], object)
										triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                                                map_triplets.Add(enhid_uriref, seckey2uri[secondarykey], object)
										map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
										if map_triplets[enhid_uriref][seckey2uri[secondarykey]][object] == 1 {
										//Because in the mapping process, different strings can use the same ID,
										//which would generate duplicate triplets
						                                        output_file.WriteString(triplet)
											set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
        	                        	               				num_ln = num_ln + 2
										}
										if ontology == "eco" {
											biolinkobject := rdf.CompU(rdf.Nss["biolink"], "EvidenceType")
											biolinktriplet := rdf.FormT(object, uris_cisreg["biolinkcat"],biolinkobject)
                                                                                	map_triplets.Add(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                                                	if map_triplets[object][uris_cisreg["biolinkcat"]][biolinkobject] == 1 {
                                                                                		output_file.WriteString(biolinktriplet)
                                                                                        	num_ln++
                                                                                	}
										}
                                                                        }
                                                                }
                                                        }
                                                }
						continue
                                        } else if secondarykey == "target_gene" {
						var gsymbgw []string
						if _, ok := genes_xmap.Lblg[tertiarykey]["bgwg"]; ok {
							gsymbgw = genes_xmap.Lblg[tertiarykey]["bgwg"].Keys()
						} else if _, ok := genes_xmap.Syng[tertiarykey]["bgwg"]; ok {
							gsymbgw = genes_xmap.Syng[tertiarykey]["bgwg"].Keys()
						} else if _, ok := genes_xmap.Ncbig[tertiarykey]["bgwg"]; ok {
							gsymbgw = genes_xmap.Ncbig[tertiarykey]["bgwg"].Keys()
						} else {
							//log.Warn("Warning: Biogateway does not contain the gene ", tertiarykey)
                                                        //fmt.Println("Warning: Biogateway does not contain the gene ", tertiarykey)
							continue}

						for _,nwgsymbgw := range gsymbgw {
                                                	object = rdf.CompU("http://rdf.biogateway.eu/gene/", nwgsymbgw)
                                                        triplet := rdf.FormT(enhid_uriref, seckey2uri[secondarykey], object)
                                                        triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                        map_triplets.Add(enhid_uriref, seckey2uri[secondarykey], object)
                                                        map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                        if map_triplets[enhid_uriref][seckey2uri[secondarykey]][object] == 1 {
                                                        	output_file.WriteString(triplet)
                                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                                num_ln = num_ln + 2
                                                	}
							if _,ok := symbol2enrich[tertiarykey]; ok {
								old_object := object
								ncbigene := symbol2enrich[tertiarykey]["ncbigene"].Keys()[0]
								//Because each gene has only one ncbigene ID
								if ncbigene != "-" {
									object = rdf.CompU("http://identifiers.org/ncbigene/", ncbigene)
									triplet = rdf.FormT(old_object, uris_cisreg["sth2exm"], object)
									map_triplets.Add(old_object, uris_cisreg["sth2exm"], object)
                                        	                	if map_triplets[old_object][uris_cisreg["sth2exm"]][object] == 1 {
                                                	                	output_file.WriteString(triplet)
                                                        	        	num_ln++
                                                        		}
									biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Gene")
        	                                          		biolinktriplet := rdf.FormT(object, uris_cisreg["biolinkcat"], biolinkobject)
                	                                        	map_triplets.Add(object, uris_cisreg["biolinkcat"], biolinkobject)
                        	                                	if map_triplets[object][uris_cisreg["biolinkcat"]][biolinkobject] == 1 {
                                	                        		output_file.WriteString(biolinktriplet)
                                        	                        	num_ln++
									}
                                                		}
							}
						}
						continue
                                        } else if secondarykey == "tfactor" {
                                                if _, ok := symbol2enrich[tertiarykey]; ok {
                                                        uniprot_ids := symbol2enrich[tertiarykey]["uniprot"].Keys()
                                                        for _,uniprot_id := range uniprot_ids {
								if upacbgw, ok := genes_xmap.Upac[uniprot_id]["bgwp"]; ok && uniprot_id != "-" {
									for _,nwupacbgw := range upacbgw.Keys() {
                                                                        	object = rdf.CompU("http://rdf.biogateway.eu/prot/", nwupacbgw)
										triplet := rdf.FormT(enhid_uriref, seckey2uri[secondarykey], object)
										triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                				map_triplets.Add(enhid_uriref, seckey2uri[secondarykey], object)
										map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
										if map_triplets[enhid_uriref][seckey2uri[secondarykey]][object] == 1 {
			                                                		output_file.WriteString(triplet)
											set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                					num_ln = num_ln + 2
										}
										triplet = rdf.FormT(object, seckey2uri[secondarykey], enhid_uriref)
										triplet_inst = rdf.FormT(object, seckey2uri[secondarykey], instid_uriref)
                                                                        	map_triplets.Add(object, seckey2uri[secondarykey], enhid_uriref)
										map_triplets.Add(object, seckey2uri[secondarykey], instid_uriref)
										if map_triplets[object][seckey2uri[secondarykey]][enhid_uriref] == 1 {
                                                                        		output_file.WriteString(triplet)
											set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                                        		num_ln = num_ln + 2

											triplet = rdf.FormT(object, uris_cisreg["sub2cls"],
													uris_cisreg["tfactor"])
											map_triplets.Add(object, uris_cisreg["sub2cls"],
													uris_cisreg["tfactor"])
	                                                                        	output_file.WriteString(triplet)
        	                                                                	num_ln++
										}
										old_object := object
										object = rdf.CompU("http://identifiers.org/uniprot/", uniprot_id)
										triplet = rdf.FormT(old_object, uris_cisreg["sth2exm"], object)
                                                                        	map_triplets.Add(old_object, uris_cisreg["sth2exm"], object)
                                                                        	if map_triplets[old_object][uris_cisreg["sth2exm"]][object] == 1 {
                                                                                	output_file.WriteString(triplet)
                                                                                	num_ln++
                                                                        	}
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Protein")
                                                        			biolinktriplet := rdf.FormT(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                        			map_triplets.Add(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                        			if map_triplets[object][uris_cisreg["biolinkcat"]][biolinkobject] == 1 {
                                                                			output_file.WriteString(biolinktriplet)
                                                                			num_ln++
                                                        			}
									}
                                                                } else if uniprot_id != "-" {
									//log.Warn("Warning: Biogateway does not contain the tfactor ", uniprot_id)
									//fmt.Println("Warning: Biogateway does not contain the gene ", tertiarykey)
								}
                                                        }
                                                }
						continue
                                        } else if secondarykey == "biosample_name" {
						if _, ok := biosamples[tertiarykey]; ok {
							biosamples_onto := biosamples[tertiarykey].Keys()
                                                        for _, onto := range biosamples_onto {
								biosample_id := biosamples[tertiarykey][onto].Keys()[0]
								if biosample_id != "-" {
									var property string
                                                                        switch{
                                                                        case onto == "CL":
                                                                                property = uris_cisreg["seq2cell"]
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Cell")
                                                                        case onto == "CLO":
                                                                                property = uris_cisreg["seq2cell_line"]
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "CellLine")
                                                                        case onto == "UBERON":
                                                                                property = uris_cisreg["seq2anato"]
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "GrossAnatomicalStructure")
									case onto == "BTO":
										switch{
										case biosamples[tertiarykey]["type"].Keys()[0] == "CL":
											property = uris_cisreg["seq2cell_line"]
										case biosamples[tertiarykey]["type"].Keys()[0] == "CT":
											property = uris_cisreg["seq2cell"]
										case biosamples[tertiarykey]["type"].Keys()[0] == "A":
											property = uris_cisreg["seq2anato"]
										}
										object = rdf.CompU("http://purl.obolibrary.org/obo/", biosample_id)
                                                                        	triplet := rdf.FormT(enhid_uriref, property, object)
                                                                        	triplet_inst := rdf.FormT(instid_uriref, property, object)
                                                                        	map_triplets.Add(enhid_uriref, property, object)
                                                                        	map_triplets.Add(instid_uriref, property, object)
                                                                        	if map_triplets[enhid_uriref][property][object] == 1 {
                                                                                	output_file.WriteString(triplet)
                                                                                	set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                                                	num_ln = num_ln + 2
                                                                        	}
										continue
									case onto == "type":
										continue
                                                                        }
									object = rdf.CompU("http://purl.obolibrary.org/obo/", biosample_id)
									triplet := rdf.FormT(enhid_uriref, property, object)
        	                                                        triplet_inst := rdf.FormT(instid_uriref, property, object)
									map_triplets.Add(enhid_uriref, property, object)
                        	                                        map_triplets.Add(instid_uriref, property, object)
									if map_triplets[enhid_uriref][property][object] == 1 {
										output_file.WriteString(triplet)
                                        	                        	set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                	                	num_ln = num_ln + 2
									}
									biolinktriplet := rdf.FormT(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                                        map_triplets.Add(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                                        if map_triplets[object][uris_cisreg["biolinkcat"]][biolinkobject] == 1 {
                                                                        	output_file.WriteString(biolinktriplet)
                                                                                num_ln++
                                                                        }
								}
							}
						}
						continue
					} else if secondarykey == "assembly" {
						object = assembly2uri[tertiarykey][0]
					} else if secondarykey == "chr" {
						object = chromosomes[tertiarykey][0]
					} else if secondarykey == "source" {
						object = sources2uri[tertiarykey][0]
					} else if secondarykey == "start" || secondarykey == "end"{
						//object = tertiarykey
						object = strings.Join([]string{rdf.FormL(tertiarykey),"^^<http://www.w3.org/2001/XMLSchema#integer>"},"")
					} else if secondarykey == "label" {
						continue
					} else if secondarykey == "crossref" {
						object = rdf.FormU(tertiarykey)
					} else {
						object = rdf.FormL(tertiarykey)
					}

					if object != "-" {
						triplet := rdf.FormT(enhid_uriref, seckey2uri[secondarykey], object)
                                        	map_triplets.Add(enhid_uriref, seckey2uri[secondarykey], object)
                                        	output_file.WriteString(triplet)

						triplet = rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
						map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
        	                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
						num_ln = num_ln + 2

						if biolinkobject != "-" {
							biolinktriplet := rdf.FormT(object, uris_cisreg["biolinkcat"], biolinkobject)
							map_triplets.Add(object, uris_cisreg["biolinkcat"], biolinkobject)
							if map_triplets[object][uris_cisreg["biolinkcat"]][biolinkobject] == 1 {
								output_file.WriteString(biolinktriplet)
								num_ln++
							}
						}
					}
                                }
                        }
                }
		output_file.WriteString(set_inst)
        }

	//Write and return
        rdfFile2export.Write([]byte(output_file.String()))
        output_file.Reset()

	return map_triplets, old2new_IDs, new2old_IDs, num_ln, nil
}

//Func to export gene coordinates to an RDF file (use the set3D as intermediary)
        //path_source : path of the table file with the information of the genes (input)
        //path_rdf : path where we want to save the RDF file (input)
	//txid : ID of the taxon (input). Example: "9606"
	//p_genes_xmap : pointer of bgw genes map from json file (input)
	//p_chr_map : pointer of chromosomes map (input)
        //map_triplets : map3D with the triplets generated for the RDF file (output)
        //num_ln : number of triplets of the RDF file (output)
func ExportGeneCoord2rdf(path_source, path_rdf, txid string, p_genes_xmap *bgw.Xmap, p_chr_map *util.Set4D) (map_triplets util.Set3D, num_ln int, err error){
        //File creation
        rdfFile2export, err := os.Create(path_rdf)
        if err != nil {
                panic(err)
        }
        defer rdfFile2export.Close()
        var output_file strings.Builder
        keys4genecoord := make(util.SliceSet)

        //Object properties
        keys4genecoord["Opys"] = []string{
                "seq2chr",
                "ins2cls",
		"sub2cls",
		"seq2version",
		"biolinkcat",
		"sth2exm",
        }

        //Annotation properties
        keys4genecoord["Apys"] = []string{
                "sth2lbl",
        }

	//Parental classes
        keys4genecoord["Prns"] = []string{
                "assembly",     //genome assembly
                "chr",          //chromosomes
        }

        //URIS
        uris := rdf.FmtURIs(keys4genecoord)

        //BGW properties to Biolink properties
        proptypes := keys4genecoord.Keys()
        for _,proptype := range proptypes {
                if proptype == "Apys" || proptype == "Opys" {
                        properties := keys4genecoord[proptype]
                        for _,property := range properties {
                                if _, ok := propbgw2propbiolink[property]; ok {
                                        triplet := rdf.FormT(uris[property], propbgw2propbiolink[property][0], propbgw2propbiolink[property][1])
                                        output_file.WriteString(triplet)
                                        num_ln++
                                }
                        }
                }
        }

        //Header of the RDF file
        num_ln = 0
        header, num_ln_header := rdf.Capita(keys4genecoord)
        output_file.WriteString(header)
        num_ln += num_ln_header

	//Datatype properties (the function rdf.Capita does not include this type of properties)
        datatype_prop := []string{
                "seq2start",
	        "seq2end",
		"seq2strand",
        }
        for _,value := range datatype_prop {
                triplet := rdf.FormT(Dpys[value][0], uris["ins2cls"], rdf.CompU(rdf.Nss["owl"], "DatatypeProperty"))
                output_file.WriteString(triplet)

                triplet = rdf.FormT(Dpys[value][0], uris["sth2lbl"], rdf.FormL(Dpys[value][1]))
                output_file.WriteString(triplet)

                triplet = rdf.FormT(Dpys[value][0], propbgw2propbiolink[value][0], propbgw2propbiolink[value][1])
                output_file.WriteString(triplet)
                num_ln = num_ln + 3
        }

        //prefLabels of chr
	chr_map := *p_chr_map
	chr_keys := chr_map[txid].Keys()
	for _,chr := range chr_keys{
		set2D := chr_map[txid][chr]
		var subject, object string

		chr_urirefseq_keys := set2D["chr_urirefseq"].Keys()
                for _,chr_urirefseq := range chr_urirefseq_keys{
                	subject = chr_urirefseq
                }
		chr_label_keys := set2D["chr_label"].Keys()
                for _,chr_label := range chr_label_keys {
                	object = chr_label
                }
		chr := strings.Join([]string{"<https://www.ncbi.nlm.nih.gov/nuccore/", subject, ">"}, "")
		triplet := rdf.FormT(chr, uris["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
		output_file.WriteString(triplet)

		triplet = rdf.FormT(chr, uris["sub2cls"], uris["chr"])
		output_file.WriteString(triplet)

                triplet = rdf.FormT(chr, uris["sth2lbl"], rdf.FormL(object))
		output_file.WriteString(triplet)

		triplet = rdf.FormT(chr, uris["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "GenomicSequenceLocalization"))
                output_file.WriteString(triplet)
                num_ln = num_ln + 4
        }

        //Map secondary keys of source to URI properties
        seckey2uri := map[string]string {
                "chr":          uris["seq2chr"],
                "start":        Dpys["seq2start"][0],
                "end":          Dpys["seq2end"][0],
		"assembly":	uris["seq2version"],
		"strand":	Dpys["seq2strand"][0],
        }

	//Map for gene enrich
        name_fields := []string{
                "symbol",
                "ncbigene",
                "uniprot",
        }
	var skip []int
        symbol2enrich, err := EnrichFile2map3D("/home/jmulero/data/enrich/symbol2enrich.tsv", name_fields, skip)
        if err != nil {
                panic(err)
        }


	//Gene coordinates to triplets
        name_fields = []string{
		"ID",
                "symbol",
		"geneid",
                "chr",
                "start",
                "end",
                "strand",
		"assembly",
        }
        skip = []int{0}
        gene_coord, err := EnrichFile2map3D(path_source, name_fields, skip)
        if err != nil {
                panic(err)
        }

        map_triplets = make(util.Set3D)
        genes_xmap := *p_genes_xmap
	gene_keys := gene_coord.Keys()
        for _,gene := range gene_keys {
		set2D := gene_coord[gene]
		gene_symbol := set2D["symbol"].Keys()[0]
		geneid := set2D["geneid"].Keys()[0]
		var gsymbgw []string
                if _, ok := genes_xmap.Lblg[gene_symbol]["bgwg"]; ok {
                	gsymbgw = genes_xmap.Lblg[gene_symbol]["bgwg"].Keys()
                } else if _, ok := genes_xmap.Syng[gene_symbol]["bgwg"]; ok {
                        gsymbgw = genes_xmap.Syng[gene_symbol]["bgwg"].Keys()
                } else if _, ok := genes_xmap.Ncbig[geneid]["bgwg"]; ok {
                        gsymbgw = genes_xmap.Ncbig[geneid]["bgwg"].Keys()
                } else {
                       //log.Warn("Warning: Biogateway does not contain the gene ", gene_symbol)
                       //fmt.Println("Warning: Biogateway does not contain the gene ", gene_symbol)
                       continue}

		for _,nwgsymbgw := range gsymbgw {
			var subject, object, predicate, triplet string
			subject = rdf.CompU("http://rdf.biogateway.eu/gene/", nwgsymbgw)
			if _,ok := symbol2enrich[gene_symbol]; ok {
                        	ncbigene := symbol2enrich[gene_symbol]["ncbigene"].Keys()[0]
				if ncbigene != "-" {
                                	object = rdf.CompU("http://identifiers.org/ncbigene/", ncbigene)
                                        triplet = rdf.FormT(subject, uris["sth2exm"], object)
                                        map_triplets.Add(object, uris["sth2exm"], object)
                                        if map_triplets[object][uris["sth2exm"]][object] == 1 {
                                        	output_file.WriteString(triplet)
                                                num_ln++
                                        }
                                        biolinkobject := rdf.CompU(rdf.Nss["biolink"], "Gene")
                                        biolinktriplet := rdf.FormT(object, uris["biolinkcat"], biolinkobject)
                                        map_triplets.Add(object, uris["biolinkcat"], biolinkobject)
                                        if map_triplets[object][uris["biolinkcat"]][biolinkobject] == 1 {
                                        	output_file.WriteString(biolinktriplet)
                                                num_ln++
                                        }
                               }
			}
			secondarykeys := set2D.Keys()
                        for _,secondarykey := range secondarykeys {
				tertiarykeys := set2D[secondarykey].Keys()
                        	for _,tertiarykey := range tertiarykeys {
                                        if secondarykey == "start" || secondarykey == "end" {
                                                //object = tertiarykey
                                               	object = strings.Join([]string{rdf.FormL(tertiarykey),"^^<http://www.w3.org/2001/XMLSchema#integer>"},"")
                                        } else if secondarykey == "chr" {
						if _, ok := chr_map[txid][tertiarykey]["chr_urirefseq"]; ok {
	                                        	for key,_ := range chr_map[txid][tertiarykey]["chr_urirefseq"]{
        	                                        	object = strings.Join([]string{"<https://www.ncbi.nlm.nih.gov/nuccore/", key, ">"}, "")
                	                                }
						} else {object = strings.Join([]string{"<https://www.ncbi.nlm.nih.gov/nuccore/", tertiarykey, ">"}, "")}
                                        } else if secondarykey == "assembly" {
						object = rdf.CompU("https://www.ncbi.nlm.nih.gov/assembly/", tertiarykey)
						triplet = rdf.FormT(object, uris["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "Genome"))
						map_triplets.Add(object, uris["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "Genome"))
						if map_triplets[object][uris["biolinkcat"]][rdf.CompU(rdf.Nss["biolink"], "Genome")] == 1 {
			                                output_file.WriteString(triplet)
                        			        num_ln++
                        			}
					} else if secondarykey == "strand" {
						object = rdf.FormL(tertiarykey)
					} else {continue}
                                        predicate = seckey2uri[secondarykey]
                                        triplet = rdf.FormT(subject, predicate, object)
                                        map_triplets.Add(subject, predicate, object)
                                        output_file.WriteString(triplet)
                                        num_ln++
                                }
                        }
		}
	}

	//write and return
        rdfFile2export.Write([]byte(output_file.String()))
        output_file.Reset()

        return map_triplets, num_ln, nil
}

//Func to export the table of TAD to RDF file (use the set3D as intermediary)
        //path_source : path of the file, the table with the TAD (input)
        //path_rdf : path where we want to save the RDF file (input)
        //name_source : name of the database used as source
        //map_triplets : map3D with the triplets generated for the RDF file
        //old2new_IDs: map3D with the old IDs and the new IDs
        //new2old_IDs: map3D with the new and the old IDs
        //num_ln : number of triplets of the RDF file
func ExportTAD2rdf(path_source, path_rdf, source_name string, p_old2new_IDs, p_new2old_IDs *util.Set3D) (map_triplets, old2new_IDs, new2old_IDs util.Set3D, num_ln int, err error){
        //File creation
        rdfFile2export, err := os.Create(path_rdf)
        if err != nil {
                panic(err)
        }
        defer rdfFile2export.Close()
        var output_file strings.Builder
        keys4tad := make(util.SliceSet)

        //Annotation properties
        keys4tad["Apys"] = []string{
                "sth2dfn",
                "sth2lbl",
                "sth2syn",
        }

        //Object properties
        keys4tad["Opys"] = []string{
                "ins2cls",
                "sub2cls",
                "sth2ori",
                "sth2evd",
                "sth2mtd",
                "sth2src",
                "gn2txn",
                "seq2chr",
                "seq2version",
		"seq2cell_line",
		"seq2cell",
		"seq2anato",
		"biolinkcat",
		"is_a",
		"sth2exm",
        }

        //Parental classes
        keys4tad["Prns"] = []string{
                "tad",
		"assembly",
		"chr",
		"db",
        }

        //Header of the RDF file
        num_ln = 0
        header, num_ln_header := rdf.Capita(keys4tad)
        output_file.WriteString(header)
        num_ln += num_ln_header

        //URIs for TAD
        uris_tad := rdf.FmtURIs(keys4tad)

	//BGW properties to Biolink properties
        proptypes := keys4tad.Keys()
        for _,proptype := range proptypes {
                if proptype == "Apys" || proptype == "Opys" {
                        properties := keys4tad[proptype]
                        for _,property := range properties {
                                if _, ok := propbgw2propbiolink[property]; ok {
                                        triplet := rdf.FormT(uris_tad[property], propbgw2propbiolink[property][0], propbgw2propbiolink[property][1])
                                        output_file.WriteString(triplet)
                                        num_ln++
                                }
                        }
                }
        }

        //Datatype properties
        datatype_prop := []string{
                "seq2start",
                "seq2end",
        }
        for _,value := range datatype_prop {
                triplet := rdf.FormT(Dpys[value][0], uris_tad["ins2cls"], rdf.CompU(rdf.Nss["owl"], "DatatypeProperty"))
                output_file.WriteString(triplet)

                triplet = rdf.FormT(Dpys[value][0], uris_tad["sth2lbl"], rdf.FormL(Dpys[value][1]))
                output_file.WriteString(triplet)

                triplet = rdf.FormT(Dpys[value][0], propbgw2propbiolink[value][0], propbgw2propbiolink[value][1])
                output_file.WriteString(triplet)
                num_ln = num_ln + 3
        }

        //AltLabel TAD and biolink category
        triplet := rdf.FormT(uris_tad["tad"], uris_tad["sth2syn"], rdf.FormL("topologically associated domain"))
        output_file.WriteString(triplet)

        triplet = rdf.FormT(uris_tad["tad"], uris_tad["is_a"], rdf.CompU(rdf.Nss["biolink"], "NucleicAcidEntity"))
        output_file.WriteString(triplet)
        num_ln = num_ln + 2

        //Chromosomes
        keys_dic := chromosomes.Keys()
        for _,value := range keys_dic {
                triplet := rdf.FormT(chromosomes[value][0], uris_tad["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
                output_file.WriteString(triplet)

                triplet = rdf.FormT(chromosomes[value][0], uris_tad["sub2cls"], chromosomes[value][2])
                output_file.WriteString(triplet)

                triplet = rdf.FormT(chromosomes[value][0], uris_tad["sth2lbl"], rdf.FormL(chromosomes[value][1]))
                output_file.WriteString(triplet)

                triplet = rdf.FormT(chromosomes[value][0], uris_tad["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "GenomicSequenceLocalization"))
                output_file.WriteString(triplet)
                num_ln = num_ln + 4
        }

	//Assemblies
        keys_dic = assembly2uri.Keys()
        for _,value := range keys_dic {
                triplet := rdf.FormT(assembly2uri[value][0], uris_tad["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
                output_file.WriteString(triplet)

                triplet = rdf.FormT(assembly2uri[value][0], uris_tad["sub2cls"], assembly2uri[value][2])
                output_file.WriteString(triplet)

                triplet = rdf.FormT(assembly2uri[value][0], uris_tad["sth2lbl"], rdf.FormL(assembly2uri[value][1]))
                output_file.WriteString(triplet)

                triplet = rdf.FormT(assembly2uri[value][0], uris_tad["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "Genome"))
                output_file.WriteString(triplet)
                num_ln = num_ln + 4
        }

        //Database
        triplet = rdf.FormT(sources2uri[source_name][0], uris_tad["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
        output_file.WriteString(triplet)

        triplet = rdf.FormT(sources2uri[source_name][0], uris_tad["sub2cls"], sources2uri[source_name][2])
        output_file.WriteString(triplet)

        triplet = rdf.FormT(sources2uri[source_name][0], uris_tad["sth2lbl"], rdf.FormL(sources2uri[source_name][1]))
        output_file.WriteString(triplet)

        triplet = rdf.FormT(sources2uri[source_name][0], uris_tad["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "InformationResource"))
        output_file.WriteString(triplet)
        num_ln = num_ln + 4


        //Source to RDF
        resource, old2new_IDs, new2old_IDs, err := TAD2map3D(path_source, p_old2new_IDs, p_new2old_IDs)
        if err != nil {
                panic(err)
	}

	//Map secondary keys of source to URI properties
        seckey2uri := map[string]string {
                "chr":          uris_tad["seq2chr"],
                "start":        Dpys["seq2start"][0],
                "end":          Dpys["seq2end"][0],
                "assembly":     uris_tad["seq2version"],
                "label":        uris_tad["sth2lbl"],
                //"orig_id":      uris_tad["sth2syn"],
                "crossref":     uris_tad["sth2ori"],
                "pubmed":       uris_tad["sth2evd"],
                "method":       uris_tad["sth2mtd"],
                "source":       uris_tad["sth2src"],
        }

        //Map 3D of methods
        name_fields := []string{
                "method_label",
                "efo",
                "obi",
                "eco",
                "bao",
                "ncit",
                "mi",
        }
        var skip []int
        methods2ids,err := EnrichFile2map3D("/home/jmulero/data/enrich/methods.tsv", name_fields, skip)
        if err != nil {
                panic(err)
        }

	//Map of method to URI properties
        methkey2uri := map[string]string{
                "efo":  rdf.Nss["efo"],
                "obi":  rdf.Nss["obo"],
                "eco":  rdf.Nss["obo"],
                "bao":  rdf.Nss["bao"],
                "ncit": rdf.Nss["obo"],
                "mi":   rdf.Nss["obo"],
        }

        //Map for biosamples
        name_fields = []string{
                "label",
                "CLO",
                "CL",
                "UBERON",
		"BTO",
		"type",
        }
        biosamples, err := EnrichFile2map3D("/home/jmulero/data/enrich/biosamples.tsv", name_fields, skip)
        if err != nil {
                panic(err)
        }

        //Set3D of the source to triplets
        map_triplets = make(util.Set3D)
        str_tad9606 := "http://rdf.biogateway.eu/tad/9606/"
	uriref_cls := rdf.CompU(rdf.Nss["owl"], "Class")

        tad_ids := resource.Keys()
        for _,primarykey := range tad_ids {
		set2D := resource[primarykey]

		//The tad is a class
                tadid_uriref := rdf.CompU(str_tad9606, primarykey)
                triplet := rdf.FormT(tadid_uriref, uris_tad["ins2cls"], uriref_cls)
                map_triplets.Add(tadid_uriref, uris_tad["ins2cls"], uriref_cls)
                output_file.WriteString(triplet)
                num_ln++

                //The tad is a subclass of TAD (class of Sequence Ontology)
                triplet = rdf.FormT(tadid_uriref, uris_tad["sub2cls"], uris_tad["tad"])
                map_triplets.Add(tadid_uriref, uris_tad["sub2cls"], uris_tad["tad"])
                output_file.WriteString(triplet)
                num_ln++

		//Definition of the class
                var chr, start, end string
                for label,_ := range set2D["label"]{
			label = strings.Replace(label, "\"", "", -1)
                        split := strings.Split(label, ":")
                        chr = split[0]
                        split = strings.Split(split[1], "-")
                        start = split[0]
                        end = split[1]
                }
                definition := fmt.Sprintf("Topologically associated domain located in Homo sapiens %s between %s and %s", chr, start, end)
                triplet = rdf.FormT(tadid_uriref, uris_tad["sth2dfn"], rdf.FormL(definition))
                map_triplets.Add(tadid_uriref, uris_tad["sth2dfn"], rdf.FormL(definition))
                output_file.WriteString(triplet)
                num_ln++

		//prefLabel of classes
                triplet = rdf.FormT(tadid_uriref, uris_tad["sth2lbl"], rdf.FormL(strings.Join([]string{"tad/", primarykey}, "")))
                map_triplets.Add(tadid_uriref, uris_tad["sth2lbl"], rdf.FormL(strings.Join([]string{"tad/", primarykey}, "")))
                output_file.WriteString(triplet)
                num_ln++

                //Instances
                var instid_uriref, database string
                for source,_ := range set2D["source"]{
			database = strings.ToLower(source)
                        instid_uriref = fmt.Sprintf("<%s%s#%s>", str_tad9606, primarykey, database)
                }
                map_triplets.Add(instid_uriref, uris_tad["ins2cls"], tadid_uriref)
                inst2class := rdf.FormT(instid_uriref, uris_tad["ins2cls"], tadid_uriref)
                set_inst := fmt.Sprintf("%s", inst2class)
                num_ln++

                //prefLabel of instances
                triplet = rdf.FormT(instid_uriref, uris_tad["sth2lbl"], rdf.FormL(strings.Join([]string{"tad#", primarykey, "#", database}, "")))
                map_triplets.Add(instid_uriref, uris_tad["sth2lbl"], rdf.FormL(strings.Join([]string{"tad#", primarykey, "#", database}, "")))
                set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

                //The sequences belong to Homo sapiens
                humantaxon_uri := rdf.CompU(rdf.Nss["obo"], "NCBITaxon_9606")
                triplet = rdf.FormT(tadid_uriref, uris_tad["gn2txn"], humantaxon_uri)
                map_triplets.Add(tadid_uriref, uris_tad["gn2txn"], humantaxon_uri)
                output_file.WriteString(triplet)

		triplet = rdf.FormT(instid_uriref, uris_tad["gn2txn"], humantaxon_uri)
                map_triplets.Add(instid_uriref, uris_tad["gn2txn"], humantaxon_uri)
                set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln = num_ln + 2

		//NCBITaxon is a category of biolink:OrganismTaxon
                biolink_taxon := rdf.CompU(rdf.Nss["biolink"], "OrganismTaxon")
                triplet = rdf.FormT(humantaxon_uri, uris_tad["biolinkcat"], biolink_taxon)
                map_triplets.Add(humantaxon_uri, uris_tad["biolinkcat"], biolink_taxon)
                if map_triplets[humantaxon_uri][uris_tad["biolinkcat"]][biolink_taxon] == 1 {
                        output_file.WriteString(triplet)
                        num_ln++
                }

                //Loop for secondary and tertiary keys
		secondarykeys := set2D.Keys()
                for _,secondarykey := range secondarykeys {
                        set1D := set2D[secondarykey]
                        tertiarykeys := set1D.Keys()
                        for _,tertiarykey := range tertiarykeys {
                                object, biolinkobject := "-", "-"
                                if  tertiarykey != "-" {
                                        if secondarykey == "pubmed" {
						object = rdf.CompU(rdf.Nss["pubmed"], tertiarykey)
                                                triplet := rdf.FormT(tadid_uriref, seckey2uri[secondarykey], object)
                                                triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(tadid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                output_file.WriteString(triplet)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln = num_ln + 2

                                                old_object := object
                                                object = rdf.CompU("http://www.ncbi.nlm.nih.gov/pubmed/", tertiarykey)
                                                triplet = rdf.FormT(old_object, uris_tad["sth2exm"], object)
                                                map_triplets.Add(old_object, uris_tad["sth2exm"], object)
                                                if map_triplets[old_object][uris_tad["sth2exm"]][object] == 1 {
                                                        output_file.WriteString(triplet)
                                                        num_ln++
                                                }

                                                biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Article")
                                                biolinktriplet := rdf.FormT(object, uris_tad["biolinkcat"], biolinkobject)
                                                map_triplets.Add(object, uris_tad["biolinkcat"], biolinkobject)
                                                if map_triplets[object][uris_tad["biolinkcat"]][biolinkobject] == 1 {
                                                        output_file.WriteString(biolinktriplet)
                                                        num_ln++
                                                }
                                                continue
                                        } else if secondarykey == "method" {
						if _, ok := methods2ids[tertiarykey]; ok {
                                                        ontologies := methods2ids[tertiarykey].Keys()
                                                        for _,ontology := range ontologies {
                                                                methods_ids := methods2ids[tertiarykey][ontology].Keys()
                                                                for _,method_id := range methods_ids {
                                                                        if method_id != "-" {
                                                                                object = rdf.CompU(methkey2uri[ontology], method_id)
                                                                                triplet := rdf.FormT(tadid_uriref, seckey2uri[secondarykey], object)
                                                                                triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                                                map_triplets.Add(tadid_uriref, seckey2uri[secondarykey], object)
                                                                                map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                                                if map_triplets[tadid_uriref][seckey2uri[secondarykey]][object] == 1 {
                                                                                //Because in the mapping process, different strings can use the same ID,
                                                                                //which would generate duplicate triplets
                                                                                        output_file.WriteString(triplet)
                                                                                        set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                                                        num_ln = num_ln + 2
                                                                                }
										if ontology == "eco" {
											biolinkobject := rdf.CompU(rdf.Nss["biolink"], "EvidenceType")
                                                                                	biolinktriplet := rdf.FormT(object, uris_tad["biolinkcat"], biolinkobject)
                                                                                	map_triplets.Add(object, uris_tad["biolinkcat"], biolinkobject)
                                                                                	if map_triplets[object][uris_tad["biolinkcat"]][biolinkobject] == 1 {
                                                                                		output_file.WriteString(biolinktriplet)
                                                                                        	num_ln++
                                                                                	}
										}
									}
                                                                }
                                                        }
                                                }
                                                continue
                                        } else if secondarykey == "biosample_name" {
                                                if _, ok := biosamples[tertiarykey]; ok {
                                                        biosamples_onto := biosamples[tertiarykey].Keys()
                                                        for _, onto := range biosamples_onto {
                                                                biosample_id := biosamples[tertiarykey][onto].Keys()[0]
                                                                if biosample_id != "-" {
                                                                        var property string
                                                                        switch{
                                                                        case onto == "CL":
                                                                                property = uris_tad["seq2cell"]
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Cell")
                                                                        case onto == "CLO":
                                                                                property = uris_tad["seq2cell_line"]
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "CellLine")
                                                                        case onto == "UBERON":
                                                                                property = uris_tad["seq2anato"]
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "GrossAnatomicalStructure")
									case onto == "BTO":
                                                                                switch{
                                                                                case biosamples[tertiarykey]["type"].Keys()[0] == "CL":
                                                                                        property = uris_tad["seq2cell_line"]
                                                                                case biosamples[tertiarykey]["type"].Keys()[0] == "CT":
                                                                                        property = uris_tad["seq2cell"]
                                                                                case biosamples[tertiarykey]["type"].Keys()[0] == "A":
                                                                                        property = uris_tad["seq2anato"]
                                                                                }
										object = rdf.CompU("http://purl.obolibrary.org/obo/", biosample_id)
                                                                        	triplet := rdf.FormT(tadid_uriref, property, object)
                                                                        	triplet_inst := rdf.FormT(instid_uriref, property, object)
                                                                        	map_triplets.Add(tadid_uriref, property, object)
                                                                        	map_triplets.Add(instid_uriref, property, object)
                                                                        	if map_triplets[tadid_uriref][property][object] == 1 {
                                                                                	output_file.WriteString(triplet)
                                                                                	set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                                                	num_ln = num_ln + 2
                                                                        	}
										continue
                                                                        case onto == "type":
                                                                                continue
                                                                        }
									object = rdf.CompU("http://purl.obolibrary.org/obo/", biosample_id)
                                                                        triplet := rdf.FormT(tadid_uriref, property, object)
                                                                        triplet_inst := rdf.FormT(instid_uriref, property, object)
                                                                        map_triplets.Add(tadid_uriref, property, object)
                                                                        map_triplets.Add(instid_uriref, property, object)
                                                                        if map_triplets[tadid_uriref][property][object] == 1 {
                                                                                output_file.WriteString(triplet)
                                                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                                                num_ln = num_ln + 2
                                                                        }
									biolinktriplet := rdf.FormT(object, uris_tad["biolinkcat"], biolinkobject)
                                                                        map_triplets.Add(object, uris_tad["biolinkcat"], biolinkobject)
                                                                        if map_triplets[object][uris_tad["biolinkcat"]][biolinkobject] == 1 {
                                                                                output_file.WriteString(biolinktriplet)
                                                                                num_ln++
                                                                        }
                                                                }
                                                        }
                                                }
                                                continue
					} else if secondarykey == "assembly" {
                                                object = assembly2uri[tertiarykey][0]
                                        } else if secondarykey == "chr" {
                                                object = chromosomes[tertiarykey][0]
                                        } else if secondarykey == "source" {
                                                object = sources2uri[tertiarykey][0]
                                        } else if secondarykey == "start" || secondarykey == "end" {	// || secondarykey == "label" {
                                                //object = tertiarykey
                                                object = strings.Join([]string{rdf.FormL(tertiarykey),"^^<http://www.w3.org/2001/XMLSchema#integer>"},"")
					} else if secondarykey == "label" {
						continue
					} else if secondarykey == "crossref" {
                                        	object = rdf.FormU(tertiarykey)
					} else {
                                                object = rdf.FormL(tertiarykey)
                                        }

					if object != "-" {
                                                triplet := rdf.FormT(tadid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(tadid_uriref, seckey2uri[secondarykey], object)
                                                output_file.WriteString(triplet)

                                                triplet = rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                                                num_ln = num_ln + 2

                                                if biolinkobject != "-" {
                                                        biolinktriplet := rdf.FormT(object, uris_tad["biolinkcat"], biolinkobject)
                                                        map_triplets.Add(object, uris_tad["biolinkcat"], biolinkobject)
                                                        if map_triplets[object][uris_tad["biolinkcat"]][biolinkobject] == 1 {
                                                                output_file.WriteString(biolinktriplet)
                                                                num_ln++
                                                        }
                                                }
					}
                                }
                        }
                }
                output_file.WriteString(set_inst)
        }

	//Write and return
        rdfFile2export.Write([]byte(output_file.String()))
        output_file.Reset()

        return map_triplets, old2new_IDs, new2old_IDs, num_ln, nil
}

func WriteIDFile(map_IDs util.Set3D, path_saveIDs string) (err error){
	//File creation
        file2export, err := os.Create(path_saveIDs)
        if err != nil {
                panic(err)
        }
        defer file2export.Close()
        var output_file strings.Builder

	//Map
	primarykeys := map_IDs.Keys()
	for _,primarykey := range primarykeys {
		secondarykeys := map_IDs[primarykey].Keys()
		for _,secondarykey := range secondarykeys {
			tertiarykeys := map_IDs[primarykey][secondarykey].Keys()
			for _,tertiarykey := range tertiarykeys {
				line := strings.Join([]string{primarykey, "\t", tertiarykey, "\n"}, "")
				output_file.WriteString(line)
			}
		}
	}

	//Write and return
	file2export.Write([]byte(output_file.String()))
	output_file.Reset()

	return nil
}

func main(){
	//Map old2new_IDs
        name_fields := []string{
                "old_ID",
                "new_ID",
        }
	var skip []int
	crm_old2new_IDs, err := EnrichFile2map3D("/home/jmulero/data/crm/crm_old2new_IDs.tsv", name_fields, skip)
	tad_old2new_IDs, err := EnrichFile2map3D("/home/jmulero/data/tad/tad_old2new_IDs.tsv", name_fields, skip)
        if err != nil {
                panic(err)
        }

        //Map new2old_IDs
        name_fields = []string{
                "new_ID",
                "old_ID",
        }
	crm_new2old_IDs, err := EnrichFile2map3D("/home/jmulero/data/crm/crm_new2old_IDs.tsv", name_fields, skip)
	tad_new2old_IDs, err := EnrichFile2map3D("/home/jmulero/data/tad/tad_new2old_IDs.tsv", name_fields, skip)
        if err != nil {
                panic(err)
        }

	//Set of cisreg sources
	enh_sources := map[string][]string{
		"ENdb":		{"/home/jmulero/databases/crm/ENdb.tsv", "/home/jmulero/rdfs/crm/ENdb.nt"},
		"EnDisease":	{"/home/jmulero/databases/crm/EnDisease.tsv", "/home/jmulero/rdfs/crm/EnDisease.nt"},
		"dbSUPER":	{"/home/jmulero/databases/crm/dbSUPER.tsv", "/home/jmulero/rdfs/crm/dbSUPER.nt"},
		"FANTOM5":	{"/home/jmulero/databases/crm/FANTOM5.tsv", "/home/jmulero/rdfs/crm/FANTOM5.nt"},
		"VISTA":	{"/home/jmulero/databases/crm/VISTA.tsv", "/home/jmulero/rdfs/crm/VISTA.nt"},
		"Ensembl":	{"/home/jmulero/databases/crm/Ensembl.tsv", "/home/jmulero/rdfs/crm/Ensembl.nt"},
		"RefSeq":	{"/home/jmulero/databases/crm/RefSeq.tsv", "/home/jmulero/rdfs/crm/RefSeq.nt"},
		"RAEdb":	{"/home/jmulero/databases/crm/RAEdb.tsv", "/home/jmulero/rdfs/crm/RAEdb.nt"},
		"EnhancerDB":	{"/home/jmulero/databases/crm/EnhancerDB.tsv", "/home/jmulero/rdfs/crm/EnhancerDB.nt"},
	}

	//Map of human genes of BGW
	genes_xmap := bgw.NewXmap()
        err = genes_xmap.Unmarshal("/home/jmulero/data/gene/9606.json") //build a map from a json
        util.CheckE(err)

	//cisreg tables to RDF
	for enh_source, path := range enh_sources{
		_, new_crm_old2new_IDs, new_crm_new2old_IDs, num_ln, err := ExportEnh2rdf(path[0], path[1], enh_source, &genes_xmap, &crm_old2new_IDs, &crm_new2old_IDs)
		if err != nil {
        	        panic(err)
	        }
		crm_old2new_IDs = new_crm_old2new_IDs
		crm_new2old_IDs = new_crm_new2old_IDs
		fmt.Printf("Number of triplets in %v: %v\n", enh_source, num_ln)
	}
	err = WriteIDFile(crm_old2new_IDs, "/home/jmulero/data/crm/crm_old2new_IDs.tsv")
        if err != nil {
        	panic(err)
        }
        err = WriteIDFile(crm_new2old_IDs, "/home/jmulero/data/crm/crm_new2old_IDs.tsv")
        if err != nil {
        	panic(err)
	}

	//Set of TAD sources
	tad_sources := map[string][]string{
		"3DGB":	{"/home/jmulero/databases/tad/3DGB.tsv", "/home/jmulero/rdfs/tad/3DGB.nt"},
		"TADKB":	{"/home/jmulero/databases/tad/TADKB.tsv", "/home/jmulero/rdfs/tad/TADKB.nt"},
	}

	//TAD tables to RDF
	for tad_source, path := range tad_sources{
		_, new_tad_old2new_IDs, new_tad_new2old_IDs, num_ln, err := ExportTAD2rdf(path[0], path[1], tad_source, &tad_old2new_IDs, &tad_new2old_IDs)
		if err != nil {
			panic(err)
		}
                tad_old2new_IDs = new_tad_old2new_IDs
                tad_new2old_IDs = new_tad_new2old_IDs
		fmt.Printf("Number of triplex in %v: %v\n", tad_source, num_ln)
	}
	err = WriteIDFile(tad_old2new_IDs, "/home/jmulero/data/tad/tad_old2new_IDs.tsv")
        if err != nil {
                panic(err)
        }
        err = WriteIDFile(tad_new2old_IDs, "/home/jmulero/data/tad/tad_new2old_IDs.tsv")
        if err != nil {
                panic(err)
	}

	//Generate chr map
	map_chr, err := Chr2map("/home/jmulero/data/gene/chromosomes.tsv")
	if err != nil {
        	panic(err)
        }

	//Genes coordinates
	files, err := ioutil.ReadDir("/home/jmulero/databases/gene/")
	if err != nil {
        	panic(err)
	}
    	for _, file := range files {
		txid := strings.Split(file.Name(), "_")[1]
		txid = strings.Split(txid, ".")[0]

		genes_xmap := bgw.NewXmap()
        	err := genes_xmap.Unmarshal(strings.Join([]string{"/home/jmulero/mappings/", txid, ".json"}, ""))
        	util.CheckE(err)

		_, num_ln, err := ExportGeneCoord2rdf(strings.Join([]string{"/home/jmulero/databases/gene/", file.Name()}, ""),
					strings.Join([]string{"/home/jmulero/rdfs/gene/gene_coord_", txid, ".nt"}, ""),
					txid, &genes_xmap, &map_chr)

		if err != nil {
                	panic(err)
        	}
      		fmt.Printf("Number of triplets in Gene coordinates of txid %v: %v\n", txid, num_ln)
	}
}


