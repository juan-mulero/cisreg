package main

//Required packages
import(
        "fmt"
        "bufio"
        "os"
	"io/ioutil"
        "strings"
        "github.com/vlmir/bgw3/src/semweb"
        "github.com/vlmir/bgw3/src/util"
	"github.com/vlmir/bgw3/src/bgw"
)

//Map of human chromosomes --> format: chromosome: {URI_chr, label_chromosome, parental_class}
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
        "chrX":	{"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000023.11>", "chr-X", "<http://purl.obolibrary.org/obo/SO_0000340>"},
        "chrY":	{"<https://www.ncbi.nlm.nih.gov/nuccore/NC_000024.10>", "chr-Y", "<http://purl.obolibrary.org/obo/SO_0000340>"},
	"chrM":	{"<https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1>", "chr-M", "<http://purl.obolibrary.org/obo/SO_0000340>"},
}

//Name diseases to URI --> name of ontology to ontology prefix to construct the complete URIs for each entitity
var diseases2uri = map[string]string{
	"DOID": rdf.Nss["doid"],
        "OMIM": rdf.Nss["omim"],
        "mesh": rdf.Nss["mesh"],
}

//Map of sources to URI --> format: database_name : {URI_database, label_database, parental_class}
var sources2uri = util.SliceSet{
	//Enhancers
        "ENdb": {"<http://www.licpathway.net/ENdb/>", "ENdb", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
        "EnDisease":    {"<http://health.tsinghua.edu.cn/jianglab/endisease/>", "EnDisease", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
        "FANTOM5":      {"<https://fantom.gsc.riken.jp/5/>", "FANTOM5", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
        "VISTA":        {"<https://enhancer.lbl.gov/>", "VISTA", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"dbSUPER":      {"<https://asntech.org/dbsuper/>", "dbSUPER", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"Ensembl":	{"<https://www.ensembl.org/Homo_sapiens/Info/Index>", "Ensembl", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"RefSeq":	{"<https://www.ncbi.nlm.nih.gov/refseq/>", "RefSeq", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"RAEdb":	{"<http://www.computationalbiology.cn/RAEdb/index.php>", "RAEdb", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"EnhancerDB":	{"<http://lcbb.swjtu.edu.cn/EnhancerDB/>", "EnhancerDB", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},

	//TADs
	"TADKB":	{"<http://dna.cs.miami.edu/TADKB/>", "TADKB", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"3DGB":	{"<http://3dgenome.fsm.northwestern.edu/index.html>", "3DGB", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},

	//Enhancers (ampliation)
	"ChromHMM":     {"<https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeBroadHmm&db=hg19>", "ChromHMM", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"CancerEnD":	{"<https://webs.iiitd.edu.in/raghava/cancerend/>", "CancerEnD", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"DiseaseEnhancer":	{"<http://biocc.hrbmu.edu.cn/DiseaseEnhancer/>", "DiseaseEnhancer", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"EnhFFL":	{"<http://lcbb.swjtu.edu.cn/EnhFFL/>", "EnhFFL", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"FOCS":	{"<http://acgt.cs.tau.ac.il/focs/>", "FOCS", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"GeneHancer":	{"<https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=geneHancer>", "GeneHancer", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"GenoSTAN":	{"<https://www.cmm.in.tum.de/public/paper/GenoSTAN/>", "GenoSTAN", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"HACER":	{"<http://bioinfo.vanderbilt.edu/AE/HACER>", "HACER", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"SEA":	{"<http://218.8.241.248:8080/SEA3/>", "SEA", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"TiED":	{"<http://lcbb.swjtu.edu.cn/TiED/>", "TiED", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"EnhancerAtlas":	{"<http://www.enhanceratlas.org/indexv2.php>", "EnhancerAtlas", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"Roadmap":	{"<https://egg2.wustl.edu/roadmap/>", "Roadmap", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"JEME":	{"<http://yiplab.cse.cuhk.edu.hk/jeme/>", "JEME", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"SCREEN":	{"<https://screen.encodeproject.org/>", "SCREEN", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"SEdb": {"<http://www.licpathway.net/sedb/>", "SEdb", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
	"scEnhancer":	{"<http://enhanceratlas.net/scenhancer/>", "scEnhancer", "<http://purl.obolibrary.org/obo/NCIT_C15426>"},
}

//Map of assembly versions --> format: assembly_name: {URI_assembly, label_assembly, alt_label_assembly, parental_class}
var assembly2uri = util.SliceSet{
        "hg38": {"<https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26>", "GRCh38", "hg38", "<http://purl.obolibrary.org/obo/SO_0001505>"},
	//"hg19": {"<https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13>", "GRCh37", "hg19", "<http://purl.obolibrary.org/obo/SO_0001505>"},
}

//Datatype properties, i.e., that have a literal as the object of the relation --> format: short_key_to_represent_property : {URI_property, label_property}
var Dpys = util.SliceSet{
	"seq2start":    {"<http://purl.obolibrary.org/obo/GENO_0000894>", "start position"},
	"seq2end":      {"<http://purl.obolibrary.org/obo/GENO_0000895>", "end position"},
}

//Map to link the BGW properties to Biolink properties --> format: key_property_hBGW: {URI_to_link_both_properties, URI_biolink_property}
var propbgw2propbiolink = util.SliceSet{
	"sth2lbl":	{"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/linkml/alias>"},
	"sub2cls":	{"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/linkml/is_a>"},
	"sth2dfn":	{"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/biolink/vocab/description>"},
	"subject":	{"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/biolink/vocab/subject>"},
	"predicate":      {"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/biolink/vocab/predicate>"},
	"object":      {"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/biolink/vocab/object>"},
	"sth2mtd":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/has_evidence>"},
	"sth2ori":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/download_url>"},
	"seq2version":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/genome_build>"},
	"sth2evd":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/publication_id>"},
	"sth2src": {"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/supporting_data_source>"},
	"tlp2tlp":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/directly_physically_interacts_with>"},
	"gn2phn":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/genotype_to_disease_association_object>"},
	"seq2txn":	{"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/biolink/vocab/in_taxon>"},
	"seq2chr":	{"<http://www.w3.org/2004/02/skos/core#exactMatch>", "<https://w3id.org/biolink/vocab/part_of>"},
	"seq2sample":	{"<http://www.w3.org/2004/02/skos/core#closeMatch>", "<https://w3id.org/biolink/vocab/supporting_study_context>"},
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

//Sequence element table from database to map3D
//We store the database in a map3D format (three keys: subject, predicate, object) to ensure consistency and unique values.
//This object will then be used to generate the rdf file.
//Input:
        //pathfile : path of the file to convert
	//graphtype: graph we want to generate. Possible values: crm, crm2gene, crm2phen, crm2tfac, tad
        //name_fields : name of the columns of the table that will be use like secondary keys
        //skip : slice with the number of the columns that we do not want to use (start in 0)
//Output:
        //map3D : map3D with de data obtained
        //err: nill if the execution is correct
func SeqElem2map3D(pathfile, graphtype string, name_fields []string, skip []int) (map3D util.Set3D, err error){
        //We open the file
        file, err := os.Open(pathfile)
        if err != nil{
                panic(err)
	}

        //Maps
	map3D = make(util.Set3D)
        genes_xmap := bgw.NewXmap()
        err = genes_xmap.Unmarshal("../../../data/gene/9606.json") //build a map from a json
        util.CheckE(err)

        //We read each line in the file
	scanner := bufio.NewScanner(file)
        i := 0
        for scanner.Scan() {
                i++
                line := scanner.Text()
                if line == "" {
                        panic(fmt.Sprintf("Empty line in the row number %v", i))
                }

                //We skip the first line because correspond with the names of the different columns
                if i != 1 {
			//We split the line in the differents fields
                        fields := strings.Split(line, "\t")

                        //We iterate these fields to create the map3D
                        //Format of the map3D --> primary key: new sequence ID; secondary key: name of the column; tertiary key: value of the column
			var primarykey string
			if graphtype == "crm" || graphtype == "tad" {
				primarykey = fields[0]
			} else if graphtype == "crm2gene" && fields[18] != "-" {
				gene := fields[18]
				if _, ok := genes_xmap.Lblg[gene]["bgwg"]; ok {
                                 	gene = genes_xmap.Lblg[gene]["bgwg"].Keys()[0]
                                } else if _, ok := genes_xmap.Syng[gene]["bgwg"]; ok {
                                        gene = genes_xmap.Syng[gene]["bgwg"].Keys()[0]
                                } else if _, ok := genes_xmap.Ncbig[gene]["bgwg"]; ok {
                                        gene = genes_xmap.Ncbig[gene]["bgwg"].Keys()[0]
				}
       	                        primarykey = strings.Join([]string{"bgw!", fields[0], "--hgncsymbol!", gene}, "")
			} else if graphtype == "crm2tfac" && fields[21] != "-" {
				primarykey = strings.Join([]string{"bgw!", fields[0], "--uniprot!", fields[21]}, "")
			} else if graphtype == "crm2phen" && fields[24] != "-" {
                                object := strings.Split(fields[24], ":")
                                primarykey = strings.Join([]string{"bgw!", fields[0], "--", strings.ToLower(object[0]), "!", object[1]}, "")
			} else {continue}

                        for j,v := range name_fields{
                                if !ExistInSlice(skip, j) {
                                        secondarykey := v
                                        tertiarykey := fields[j]
                                        map3D.Add(primarykey, secondarykey, tertiarykey)
                                }
                        }
                }
        }
        return map3D, nil
}

//Chromosomes table to map4D. Format --> primarykey: taxon_id, secondarykey: chr, tertiarykey: column_file, quaternarykey: column_value
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
			"assembly_altlabel",
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

//Func to convert one file with enrichment information to map 3D (file without header)
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

//Func to export CRM data to an RDF file (use the map3D as intermediary) --> graph crm
//Input:
        //path_source: path of the file, the table with the CRM database
        //path_rdf: path where we want to save the RDF file
        //name_source: name of the database used as source

//Output:
        //map_triplets : map3D with the CRM database data
        //num_ln : number of triplets of the RDF file
        //err: nill if the execution is correct
func ExportCRM2rdf(path_source, path_rdf, source_name string)(map_triplets util.Set3D, num_ln int, err error){
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
		"seq2txn",
		"seq2chr",
                "seq2version",
		"seq2sample",
		"biolinkcat",
		"is_a",
		"sth2clm",
		"sth2exm",
        }

        //Parental classes
        keys4cisreg["Prns"] = []string{
                "crm",		//cis-regulatory module
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

	triplet = rdf.FormT(uris_cisreg["crm"], uris_cisreg["is_a"], rdf.CompU(rdf.Nss["biolink"], "RegulatoryRegion"))
        output_file.WriteString(triplet)
        num_ln = num_ln + 2

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

                triplet = rdf.FormT(assembly2uri[value][0], uris_cisreg["sub2cls"], assembly2uri[value][3])
                output_file.WriteString(triplet)

                triplet = rdf.FormT(assembly2uri[value][0], uris_cisreg["sth2lbl"], rdf.FormL(assembly2uri[value][1]))
                output_file.WriteString(triplet)

		triplet = rdf.FormT(assembly2uri[value][0], uris_cisreg["sth2syn"], rdf.FormL(assembly2uri[value][2]))
                output_file.WriteString(triplet)

		triplet = rdf.FormT(assembly2uri[value][0], uris_cisreg["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "Genome"))
		output_file.WriteString(triplet)
                num_ln = num_ln + 5
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
	name_fields := []string{
                        "crm_ID",
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

        skip := []int{0,1,2,3,4,9,10,11,16,18,19,20,21,22,23,24,25,26,27,28,29}
        resource, err := SeqElem2map3D(path_source, "crm", name_fields, skip)
        if err != nil {
                panic(err)
        }

        //Map secondary keys of source to URI properties
        seckey2uri := map[string]string {
                "chr":		uris_cisreg["seq2chr"],
                "start":	Dpys["seq2start"][0],
                "end":		Dpys["seq2end"][0],
                "assembly":	uris_cisreg["seq2version"],
                "crossref":	uris_cisreg["sth2ori"],
                "pubmed":	uris_cisreg["sth2evd"],
                "method":	uris_cisreg["sth2mtd"],
                "source":	uris_cisreg["sth2src"],
		"biosample_name":	uris_cisreg["seq2sample"],
        }

        //Map 3D of methods
	name_fields = []string{
		"method_label",
		"efo",
		"obi",
		"eco",
		"bao",
		"ncit",
		"mi",
	}
	skip = []int{}
        methods2ids,err := EnrichFile2map3D("../../../data/enrich/methods.tsv", name_fields, skip)
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
	biosamples, err := EnrichFile2map3D("../../../data/enrich/biosamples.tsv", name_fields, skip)
	if err != nil {
                panic(err)
        }

        //Set3D of the source to triplets
	map_triplets = make(util.Set3D)
	str_cisreg9606 := "http://rdf.biogateway.eu/crm/9606/"
	uriref_cls := rdf.CompU(rdf.Nss["owl"], "Class")
	crm_ids := resource.Keys()
        for _,primarykey := range crm_ids {
		set2D := resource[primarykey]

                //The CRM is a class
                crmid_uriref := rdf.CompU(str_cisreg9606, primarykey)
                triplet := rdf.FormT(crmid_uriref, uris_cisreg["ins2cls"], uriref_cls)
		map_triplets.Add(crmid_uriref, uris_cisreg["ins2cls"], uriref_cls)
                output_file.WriteString(triplet)
                num_ln++

                //The CRM is a subclass of CRM
                triplet = rdf.FormT(crmid_uriref, uris_cisreg["sub2cls"], uris_cisreg["crm"])
		map_triplets.Add(crmid_uriref, uris_cisreg["sub2cls"], uris_cisreg["crm"])
                output_file.WriteString(triplet)
		num_ln++

		//Definition of the class
		chr := set2D["chr"].Keys()[0]
		start := set2D["start"].Keys()[0]
		end := set2D["end"].Keys()[0]
		definition := fmt.Sprintf("Cis-regulatory module located in Homo sapiens %s between %s and %s", chr, start, end)
		triplet = rdf.FormT(crmid_uriref, uris_cisreg["sth2dfn"], rdf.FormL(definition))
		map_triplets.Add(crmid_uriref, uris_cisreg["sth2dfn"], rdf.FormL(definition))
                output_file.WriteString(triplet)
                num_ln++

		//prefLabel of classes
		triplet = rdf.FormT(crmid_uriref, uris_cisreg["sth2lbl"], rdf.FormL(strings.Join([]string{"crm/", primarykey}, "")))
                map_triplets.Add(crmid_uriref, uris_cisreg["sth2lbl"], rdf.FormL(strings.Join([]string{"crm/", primarykey}, "")))
                output_file.WriteString(triplet)
                num_ln++

		//Instances
		source := set2D["source"].Keys()[0]
		instid_uriref := fmt.Sprintf("<%s%s#%s>", str_cisreg9606, primarykey, strings.ToLower(source))
                map_triplets.Add(instid_uriref, uris_cisreg["ins2cls"], crmid_uriref)
                inst2class := rdf.FormT(instid_uriref, uris_cisreg["ins2cls"], crmid_uriref)
                set_inst := fmt.Sprintf("%s", inst2class)
		num_ln++

		//Definition of the instances
		definition = fmt.Sprintf("%s according to %s", definition, source)
		triplet = rdf.FormT(instid_uriref, uris_cisreg["sth2dfn"], rdf.FormL(definition))
                map_triplets.Add(instid_uriref, uris_cisreg["sth2dfn"], rdf.FormL(definition))
		set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
		num_ln++

		//prefLabel of instances
		triplet = rdf.FormT(instid_uriref, uris_cisreg["sth2lbl"], rdf.FormL(strings.Join([]string{"crm#", primarykey, "#", strings.ToLower(source)}, "")))
                map_triplets.Add(instid_uriref, uris_cisreg["sth2lbl"], rdf.FormL(strings.Join([]string{"crm#", primarykey, "#", strings.ToLower(source)}, "")))
		set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

		//The sequences belong to Homo sapiens
		humantaxon_uri := rdf.CompU(rdf.Nss["obo"], "NCBITaxon_9606")
		triplet = rdf.FormT(crmid_uriref, uris_cisreg["seq2txn"], humantaxon_uri)
	        map_triplets.Add(crmid_uriref, uris_cisreg["seq2txn"], humantaxon_uri)
                output_file.WriteString(triplet)
		num_ln++

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
						triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
						num_ln++

                                                old_object := object
                                                object = rdf.CompU("http://www.ncbi.nlm.nih.gov/pubmed/", tertiarykey)
                                                triplet = rdf.FormT(old_object, uris_cisreg["sth2exm"], object)
                                                map_triplets.Add(old_object, uris_cisreg["sth2exm"], object)
                                                if map_triplets[old_object][uris_cisreg["sth2exm"]][object] == 1 {
							set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                                                        num_ln++
                                                }

						biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Article")
						biolinktriplet := rdf.FormT(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                map_triplets.Add(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                if map_triplets[object][uris_cisreg["biolinkcat"]][biolinkobject] == 1 {
							set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
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
										triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
										map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
										if map_triplets[instid_uriref][seckey2uri[secondarykey]][object] == 1 {
										//Because in the mapping process, different strings can use the same ID,
										//which would generate duplicate triplets
											set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
        	                        	               				num_ln++
										}
										if ontology == "eco" {
											biolinkobject := rdf.CompU(rdf.Nss["biolink"], "EvidenceType")
											biolinktriplet := rdf.FormT(object, uris_cisreg["biolinkcat"],biolinkobject)
                                                                                	map_triplets.Add(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                                                	if map_triplets[object][uris_cisreg["biolinkcat"]][biolinkobject] == 1 {
												set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                                                        	num_ln++
                                                                                	}
										}
                                                                        }
                                                                }
                                                        }
                                                } else {

							fmt.Println("Warning: unmapped method ", tertiarykey)
						}
						continue
                                        } else if secondarykey == "biosample_name" {
						if _, ok := biosamples[tertiarykey]; ok {
							biosamples_onto := biosamples[tertiarykey].Keys()
                                                        for _, onto := range biosamples_onto {
								biosample_id := biosamples[tertiarykey][onto].Keys()[0]
								if biosample_id != "-" {
									property := seckey2uri[secondarykey]
                                                                        switch{
                                                                        case onto == "CL":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Cell")
                                                                        case onto == "CLO":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "CellLine")
                                                                        case onto == "UBERON":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "GrossAnatomicalStructure")
									case onto == "BTO":
										switch{
										case biosamples[tertiarykey]["type"].Keys()[0] == "CL":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "CellLine")
										case biosamples[tertiarykey]["type"].Keys()[0] == "CT":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Cell")
										case biosamples[tertiarykey]["type"].Keys()[0] == "A":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "GrossAnatomicalStructure")
										}
									case onto == "type":
										continue
                                                                        }
									object = rdf.CompU("http://purl.obolibrary.org/obo/", biosample_id)
        	                                                        triplet_inst := rdf.FormT(instid_uriref, property, object)
                        	                                        map_triplets.Add(instid_uriref, property, object)
									if map_triplets[instid_uriref][property][object] == 1 {
                                        	                        	set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                	                	num_ln++
									}
									biolinktriplet := rdf.FormT(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                                        map_triplets.Add(object, uris_cisreg["biolinkcat"], biolinkobject)
                                                                        if map_triplets[object][uris_cisreg["biolinkcat"]][biolinkobject] == 1 {
										set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                                                num_ln++
                                                                        }
								}
							}
						} else {
							fmt.Println("Warning: unmapped biosample ", tertiarykey)
						}
						continue
					} else if secondarykey == "assembly" {
						object = assembly2uri[tertiarykey][0]
					} else if secondarykey == "chr" {
						object = chromosomes[tertiarykey][0]
					} else if secondarykey == "source" {
						object = sources2uri[tertiarykey][0]
						triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
						map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
						set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
						num_ln++
						continue
					} else if secondarykey == "start" || secondarykey == "end"{
						object = strings.Join([]string{rdf.FormL(tertiarykey),"^^<http://www.w3.org/2001/XMLSchema#integer>"},"")
					} else if secondarykey == "crossref" {
						object = rdf.FormU(tertiarykey)
						triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++
                                                continue
					} else {
						object = rdf.FormL(tertiarykey)
					}

					if object != "-" {
						triplet := rdf.FormT(crmid_uriref, seckey2uri[secondarykey], object)
                                        	map_triplets.Add(crmid_uriref, seckey2uri[secondarykey], object)
                                        	output_file.WriteString(triplet)
						num_ln++

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

	return map_triplets, num_ln, nil
}

//Func to export crm-phenotype data to an RDF file (use the map3D as intermediary) --> graph crm2phen
//Input:
        //path_source: path of the file, the table with the CRM database
        //path_rdf: path where we want to save the RDF file
        //name_source: name of the database used as source

//Output:
        //map_triplets : map3D with the CRM database data
        //num_ln : number of triplets of the RDF file
        //err: nill if the execution is correct
func ExportCRM2phen2rdf(path_source, path_rdf, source_name string)(map_triplets util.Set3D, num_ln int, err error){
	//File creation
        rdfFile2export, err := os.Create(path_rdf)
        if err != nil {
                panic(err)
        }
        defer rdfFile2export.Close()
        var output_file strings.Builder
        keys4graph := make(util.SliceSet)

        //Annotation properties
        keys4graph["Apys"] = []string{
                "sth2dfn",
                "sth2lbl",
        }

        //Object properties
        keys4graph["Opys"] = []string{
                "ins2cls",
                "sub2cls",
                "subject",
		"predicate",
		"object",
		"gn2phn",
		"sth2evd",
		"sth2mtd",
		"sth2ori",
		"sth2src",
		"biolinkcat",
		"sth2exm",
        }

        //Parental classes
        keys4graph["Prns"] = []string{
		"stm",
		"db",
		"crm2phen",
        }

        //Header of the RDF file with the definition of the properties and the parental classes
        num_ln = 0
        header, num_ln_header := rdf.Capita(keys4graph)
        output_file.WriteString(header)
        num_ln += num_ln_header

	//URIs for graph
        uris_graph := rdf.FmtURIs(keys4graph)

	//BGW properties to Biolink properties
	proptypes := keys4graph.Keys()
	for _,proptype := range proptypes {
		if proptype == "Apys" || proptype == "Opys" {
			properties := keys4graph[proptype]
			for _,property := range properties {
				if _, ok := propbgw2propbiolink[property]; ok {
					triplet := rdf.FormT(uris_graph[property], propbgw2propbiolink[property][0], propbgw2propbiolink[property][1])
					output_file.WriteString(triplet)
					num_ln++
				}
			}
		}
	}

	//Database
	triplet := rdf.FormT(sources2uri[source_name][0], uris_graph["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
        output_file.WriteString(triplet)

	triplet = rdf.FormT(sources2uri[source_name][0], uris_graph["sub2cls"], sources2uri[source_name][2])
        output_file.WriteString(triplet)

	triplet = rdf.FormT(sources2uri[source_name][0], uris_graph["sth2lbl"], rdf.FormL(sources2uri[source_name][1]))
	output_file.WriteString(triplet)

        triplet = rdf.FormT(sources2uri[source_name][0], uris_graph["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "InformationResource"))
	output_file.WriteString(triplet)
        num_ln = num_ln + 4

	//source to RDF
        name_fields := []string{
                        "subject",
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
                        "object",
                        "pubmed",
                        "method",
                        "refsnp",
                        "pubmed",
                        "method",
                        }

        skip := []int{1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,18,19,20,21,22,23,27,28,29}

        resource, err := SeqElem2map3D(path_source, "crm2phen", name_fields, skip)
        if err != nil {
                panic(err)
        }

        //Map secondary keys of source to URI properties
        seckey2uri := map[string]string {
                "subject":	uris_graph["subject"],
                "object":	uris_graph["object"],
                "crossref":	uris_graph["sth2ori"],
                "pubmed":	uris_graph["sth2evd"],
                "method":	uris_graph["sth2mtd"],
                "source":	uris_graph["sth2src"],
        }

        //Map 3D of methods
	name_fields = []string{
		"method_label",
		"efo",
		"obi",
		"eco",
		"bao",
		"ncit",
		"mi",
	}
	skip = []int{}
        methods2ids,err := EnrichFile2map3D("../../../data/enrich/methods.tsv", name_fields, skip)
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

        //Set3D of the source to triplets
	map_triplets = make(util.Set3D)
	url_graph := "http://rdf.biogateway.eu/crm2phen/"
	url_crm := "http://rdf.biogateway.eu/crm/9606/"
	uriref_cls := rdf.CompU(rdf.Nss["owl"], "Class")
	graph_ids := resource.Keys()
        for _,primarykey := range graph_ids {
		set2D := resource[primarykey]
		crm := rdf.CompU(url_crm, set2D["subject"].Keys()[0])

                //The relation is a class
                uriref_relation := rdf.CompU(url_graph, primarykey)
                triplet := rdf.FormT(uriref_relation, uris_graph["ins2cls"], uriref_cls)
		map_triplets.Add(uriref_relation, uris_graph["ins2cls"], uriref_cls)
                output_file.WriteString(triplet)
                num_ln++

                //The relation is a subclass
                triplet = rdf.FormT(uriref_relation, uris_graph["sub2cls"], uris_graph["stm"])
		map_triplets.Add(uriref_relation, uris_graph["sub2cls"], uris_graph["stm"])
                output_file.WriteString(triplet)
		num_ln++

		//The relation is a category of Biolink
                triplet = rdf.FormT(uriref_relation, uris_graph["biolinkcat"], uris_graph["crm2phen"])
                map_triplets.Add(uriref_relation, uris_graph["biolinkcat"], uris_graph["crm2phen"])
                output_file.WriteString(triplet)
                num_ln++

                //The relation has a predicate
                triplet = rdf.FormT(uriref_relation, uris_graph["predicate"], uris_graph["gn2phn"])
                map_triplets.Add(uriref_relation, uris_graph["predicate"], uris_graph["gn2phn"])
                output_file.WriteString(triplet)
                num_ln++

		//prefLabel of the relation
		triplet = rdf.FormT(uriref_relation, uris_graph["sth2lbl"], rdf.FormL(primarykey))
		map_triplets.Add(uriref_relation, uris_graph["sth2lbl"], rdf.FormL(primarykey))
                output_file.WriteString(triplet)
                num_ln++

		//Definition of the relation
                subject := set2D["subject"].Keys()[0]
                object := set2D["object"].Keys()[0]
		definition := fmt.Sprintf("Association between cis-regulatory module %s and disease %s", subject, object)
		triplet = rdf.FormT(uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
		map_triplets.Add(uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
                output_file.WriteString(triplet)
                num_ln++

		//Intance of the relation class
		database := set2D["source"].Keys()[0]
		inst_id_relation := strings.Join([]string{primarykey, "#", strings.ToLower(database)}, "")
		inst_uriref_relation := rdf.CompU(url_graph, inst_id_relation)
                map_triplets.Add(inst_uriref_relation, uris_graph["ins2cls"], uriref_relation)
                inst2class := rdf.FormT(inst_uriref_relation, uris_graph["ins2cls"], uriref_relation)
                set_inst := fmt.Sprintf("%s", inst2class)
		num_ln++

		//prefLabel of instance
		triplet = rdf.FormT(inst_uriref_relation, uris_graph["sth2lbl"], rdf.FormL(inst_id_relation))
                map_triplets.Add(inst_uriref_relation, uris_graph["sth2lbl"], inst_id_relation)
		set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

                //Definition of the instance
                definition = fmt.Sprintf("%s according to %s", definition, database)
                triplet = rdf.FormT(inst_uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
                map_triplets.Add(inst_uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
		set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

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
						triplet_inst := rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++

                                                old_object := object
                                                object = rdf.CompU("http://www.ncbi.nlm.nih.gov/pubmed/", tertiarykey)
                                                triplet = rdf.FormT(old_object, uris_graph["sth2exm"], object)
                                                map_triplets.Add(old_object, uris_graph["sth2exm"], object)
                                                if map_triplets[old_object][uris_graph["sth2exm"]][object] == 1 {
							set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                                                        num_ln++
                                                }

						biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Article")
						biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"], biolinkobject)
                                                map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                                                if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
							set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                        num_ln++
                                                }
					} else if secondarykey == "method" {
                                                if _, ok := methods2ids[tertiarykey]; ok {
							ontologies := methods2ids[tertiarykey].Keys()
                                                        for _,ontology := range ontologies {
								methods_ids := methods2ids[tertiarykey][ontology].Keys()
                                                                for _,method_id := range methods_ids {
                                                                        if method_id != "-" {
										object = rdf.CompU(methkey2uri[ontology], method_id)
										triplet_inst := rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
										map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
										if map_triplets[inst_uriref_relation][seckey2uri[secondarykey]][object] == 1 {
										//Because in the mapping process, different strings can use the same ID,
										//which would generate duplicate triplets
											set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
        	                        	               				num_ln++
										}
										if ontology == "eco" {
											biolinkobject := rdf.CompU(rdf.Nss["biolink"], "EvidenceType")
											biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"],biolinkobject)
                                                                                	map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                                                                                	if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
												set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                                                        	num_ln++
                                                                                	}
										}
                                                                        }
                                                                }
                                                        }
                                                } else {
							fmt.Println("Warning: unmapped method ", tertiarykey)
						}
					} else if secondarykey == "object" {
						split := strings.Split(tertiarykey, ":")
						if uri, ok := diseases2uri[split[0]]; ok {
							object = rdf.CompU(uri, split[1])
							triplet := rdf.FormT(uriref_relation, seckey2uri[secondarykey], object)
                                                        map_triplets.Add(uriref_relation, seckey2uri[secondarykey], object)
							output_file.WriteString(triplet)
							num_ln++

							triplet = rdf.FormT(crm, uris_graph["gn2phn"], object)
							map_triplets.Add(crm, uris_graph["gn2phn"], object)
							output_file.WriteString(triplet)
                                                        num_ln++

							if split[0] == "mesh" {
								new_object := rdf.CompU("http://identifiers.org/mesh/", split[1])
								triplet = rdf.FormT(object, uris_graph["sth2exm"], new_object)
								map_triplets.Add(object, uris_graph["sth2exm"], new_object)
								if map_triplets[object][uris_graph["sth2exm"]][new_object] == 1 {
									output_file.WriteString(triplet)
                                                                	num_ln++
								}
							} else if split[0] == "OMIM" {
                                                                new_object := rdf.CompU("http://purl.obolibrary.org/obo/OMIM_", split[1])
                                                                triplet = rdf.FormT(object, uris_graph["sth2exm"], new_object)
                                                                map_triplets.Add(object, uris_graph["sth2exm"], new_object)
                                                                if map_triplets[object][uris_graph["sth2exm"]][new_object] == 1 {
                                                                        output_file.WriteString(triplet)
                                                                        num_ln++
                                                                }
								object = new_object
							}
                                                        biolinkobject := rdf.CompU(rdf.Nss["biolink"], "Disease")
                                                        biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"], biolinkobject)
                                                        map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                                                        if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
                                                                output_file.WriteString(biolinktriplet)
                                                                num_ln++
                                                        }
						}
                                        } else if secondarykey == "subject" {
						triplet := rdf.FormT(uriref_relation, seckey2uri[secondarykey], crm)
                                                map_triplets.Add(uriref_relation, seckey2uri[secondarykey], crm)
                                                output_file.WriteString(triplet)
                                                num_ln++
					} else if secondarykey == "source" {
						object = sources2uri[tertiarykey][0]
						triplet_inst := rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
						set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++
					} else if secondarykey == "crossref" {
						object = rdf.FormU(tertiarykey)
                                                triplet_inst := rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
						set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++
					}
                                }
                        }
                }
		output_file.WriteString(set_inst)
        }

	//Write and return
        rdfFile2export.Write([]byte(output_file.String()))
        output_file.Reset()

	return map_triplets, num_ln, nil
}

//Func to export CRM-TranscriptionFactor data to an RDF file (use the map3D as intermediary) --> graph crm2tfac
//Input:
        //path_source: path of the file, the table with the CRM database
        //path_rdf: path where we want to save the RDF file
        //name_source: name of the database used as source

//Output:
        //map_triplets : map3D with the CRM database data
        //num_ln : number of triplets of the RDF file
        //err: nill if the execution is correct
func ExportCRM2tfac2rdf(path_source, path_rdf, source_name string)(map_triplets util.Set3D, num_ln int, err error){
	//File creation
        rdfFile2export, err := os.Create(path_rdf)
        if err != nil {
                panic(err)
        }
        defer rdfFile2export.Close()
        var output_file strings.Builder
        keys4graph := make(util.SliceSet)

        //Annotation properties
        keys4graph["Apys"] = []string{
                "sth2dfn",
                "sth2lbl",
		"sth2syn",
        }

        //Object properties
        keys4graph["Opys"] = []string{
                "ins2cls",
                "sub2cls",
                "subject",
		"predicate",
		"object",
		"tlp2tlp",
		"sth2evd",
		"sth2mtd",
		"sth2ori",
		"sth2src",
		"biolinkcat",
		"sth2exm",
		"seq2sample",
        }

        //Parental classes
        keys4graph["Prns"] = []string{
		"stm",
		"db",
		"tfactor",
        }

        //Header of the RDF file with the definition of the properties and the parental classes
        num_ln = 0
        header, num_ln_header := rdf.Capita(keys4graph)
        output_file.WriteString(header)
        num_ln += num_ln_header

	//URIs for graph
        uris_graph := rdf.FmtURIs(keys4graph)

	//AltLabel TF
        triplet := rdf.FormT(uris_graph["tfactor"], uris_graph["sth2syn"], rdf.FormL("TF"))
        output_file.WriteString(triplet)
        num_ln++

	//BGW properties to Biolink properties
	proptypes := keys4graph.Keys()
	for _,proptype := range proptypes {
		if proptype == "Apys" || proptype == "Opys" {
			properties := keys4graph[proptype]
			for _,property := range properties {
				if _, ok := propbgw2propbiolink[property]; ok {
					triplet := rdf.FormT(uris_graph[property], propbgw2propbiolink[property][0], propbgw2propbiolink[property][1])
					output_file.WriteString(triplet)
					num_ln++
				}
			}
		}
	}

	//Database
	triplet = rdf.FormT(sources2uri[source_name][0], uris_graph["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
        output_file.WriteString(triplet)

	triplet = rdf.FormT(sources2uri[source_name][0], uris_graph["sub2cls"], sources2uri[source_name][2])
        output_file.WriteString(triplet)

	triplet = rdf.FormT(sources2uri[source_name][0], uris_graph["sth2lbl"], rdf.FormL(sources2uri[source_name][1]))
	output_file.WriteString(triplet)

        triplet = rdf.FormT(sources2uri[source_name][0], uris_graph["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "InformationResource"))
	output_file.WriteString(triplet)
        num_ln = num_ln + 4

	//source to RDF
        name_fields := []string{
                        "subject",
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
                        "object",
                        "pubmed",
                        "method",
                        "disease",
                        "pubmed",
                        "method",
                        "refsnp",
                        "pubmed",
                        "method",
                        }

        skip := []int{1,2,3,4,5,6,7,8,9,10,11,13,15,16,18,19,20,24,25,26,27,28,29}
        resource, err := SeqElem2map3D(path_source, "crm2tfac", name_fields, skip)
        if err != nil {
                panic(err)
        }

        //Map secondary keys of source to URI properties
        seckey2uri := map[string]string {
                "subject":	uris_graph["subject"],
                "object":	uris_graph["object"],
                "crossref":	uris_graph["sth2ori"],
                "pubmed":	uris_graph["sth2evd"],
                "method":	uris_graph["sth2mtd"],
                "source":	uris_graph["sth2src"],
		"biosample_name":	uris_graph["seq2sample"],
        }

        //Map 3D of methods
	name_fields = []string{
		"method_label",
		"efo",
		"obi",
		"eco",
		"bao",
		"ncit",
		"mi",
	}
	skip = []int{}
        methods2ids,err := EnrichFile2map3D("../../../data/enrich/methods.tsv", name_fields, skip)
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
        biosamples, err := EnrichFile2map3D("../../../data/enrich/biosamples.tsv", name_fields, skip)
        if err != nil {
                panic(err)
        }

        //Set3D of the source to triplets
	map_triplets = make(util.Set3D)
	url_graph := "http://rdf.biogateway.eu/crm2tfac/"
	url_crm := "http://rdf.biogateway.eu/crm/9606/"
	uriref_cls := rdf.CompU(rdf.Nss["owl"], "Class")
	graph_ids := resource.Keys()
        for _,primarykey := range graph_ids {
		set2D := resource[primarykey]
		crm := rdf.CompU(url_crm, set2D["subject"].Keys()[0])

                //The relation is a class
                uriref_relation := rdf.CompU(url_graph, primarykey)
                triplet := rdf.FormT(uriref_relation, uris_graph["ins2cls"], uriref_cls)
		map_triplets.Add(uriref_relation, uris_graph["ins2cls"], uriref_cls)
                output_file.WriteString(triplet)
                num_ln++

                //The relation is a subclass
                triplet = rdf.FormT(uriref_relation, uris_graph["sub2cls"], uris_graph["stm"])
		map_triplets.Add(uriref_relation, uris_graph["sub2cls"], uris_graph["stm"])
                output_file.WriteString(triplet)
		num_ln++

                //The relation has a predicate
                triplet = rdf.FormT(uriref_relation, uris_graph["predicate"], uris_graph["tlp2tlp"])
                map_triplets.Add(uriref_relation, uris_graph["predicate"], uris_graph["tlp2tlp"])
                output_file.WriteString(triplet)
                num_ln++

		//prefLabel of the relation
		triplet = rdf.FormT(uriref_relation, uris_graph["sth2lbl"], rdf.FormL(primarykey))
		map_triplets.Add(uriref_relation, uris_graph["sth2lbl"], rdf.FormL(primarykey))
                output_file.WriteString(triplet)
                num_ln++

		//Definition of the relation
                subject := set2D["subject"].Keys()[0]
                object := set2D["object"].Keys()[0]
		definition := fmt.Sprintf("Association between cis-regulatory module %s and transcription factor %s", subject, object)
		triplet = rdf.FormT(uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
		map_triplets.Add(uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
                output_file.WriteString(triplet)
                num_ln++

		//Intance of the relation class
		database := set2D["source"].Keys()[0]
		inst_id_relation := strings.Join([]string{primarykey, "#", strings.ToLower(database)}, "")
		inst_uriref_relation := rdf.CompU(url_graph, inst_id_relation)
                map_triplets.Add(inst_uriref_relation, uris_graph["ins2cls"], uriref_relation)
                inst2class := rdf.FormT(inst_uriref_relation, uris_graph["ins2cls"], uriref_relation)
                set_inst := fmt.Sprintf("%s", inst2class)
		num_ln++

		//prefLabel of instance
		triplet = rdf.FormT(inst_uriref_relation, uris_graph["sth2lbl"], rdf.FormL(inst_id_relation))
                map_triplets.Add(inst_uriref_relation, uris_graph["sth2lbl"], inst_id_relation)
		set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

                //Definition of the instance
                definition = fmt.Sprintf("%s according to %s", definition, database)
                triplet = rdf.FormT(inst_uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
                map_triplets.Add(inst_uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
		set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

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
						triplet_inst := rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++

                                                old_object := object
                                                object = rdf.CompU("http://www.ncbi.nlm.nih.gov/pubmed/", tertiarykey)
                                                triplet = rdf.FormT(old_object, uris_graph["sth2exm"], object)
                                                map_triplets.Add(old_object, uris_graph["sth2exm"], object)
                                                if map_triplets[old_object][uris_graph["sth2exm"]][object] == 1 {
							set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                                                        num_ln++
                                                }

						biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Article")
						biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"], biolinkobject)
                                                map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                                                if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
							set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                        num_ln++
                                                }
					} else if secondarykey == "method" {
                                                if _, ok := methods2ids[tertiarykey]; ok {
							ontologies := methods2ids[tertiarykey].Keys()
                                                        for _,ontology := range ontologies {
								methods_ids := methods2ids[tertiarykey][ontology].Keys()
                                                                for _,method_id := range methods_ids {
                                                                        if method_id != "-" {
										object = rdf.CompU(methkey2uri[ontology], method_id)
										triplet_inst := rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
										map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
										if map_triplets[inst_uriref_relation][seckey2uri[secondarykey]][object] == 1 {
										//Because in the mapping process, different strings can use the same ID,
										//which would generate duplicate triplets
											set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
        	                        	               				num_ln++
										}
										if ontology == "eco" {
											biolinkobject := rdf.CompU(rdf.Nss["biolink"], "EvidenceType")
											biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"],biolinkobject)
                                                                                	map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                                                                                	if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
												set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                                                        	num_ln++
                                                                                	}
										}
                                                                        }
                                                                }
                                                        }
                                                } else {
							fmt.Println("Warning: unmapped method ", tertiarykey)
						}
                                        } else if secondarykey == "biosample_name" {
						if _, ok := biosamples[tertiarykey]; ok {
							biosamples_onto := biosamples[tertiarykey].Keys()
                                                        for _, onto := range biosamples_onto {
								biosample_id := biosamples[tertiarykey][onto].Keys()[0]
								if biosample_id != "-" {
									property := seckey2uri[secondarykey]
                                                                        switch{
                                                                        case onto == "CL":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Cell")
                                                                        case onto == "CLO":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "CellLine")
                                                                        case onto == "UBERON":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "GrossAnatomicalStructure")
									case onto == "BTO":
										switch{
										case biosamples[tertiarykey]["type"].Keys()[0] == "CL":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "CellLine")
										case biosamples[tertiarykey]["type"].Keys()[0] == "CT":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Cell")
										case biosamples[tertiarykey]["type"].Keys()[0] == "A":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "GrossAnatomicalStructure")
										}
									case onto == "type":
										continue
                                                                        }
									object = rdf.CompU("http://purl.obolibrary.org/obo/", biosample_id)
        	                                                        triplet_inst := rdf.FormT(inst_uriref_relation, property, object)
                        	                                        map_triplets.Add(inst_uriref_relation, property, object)
									if map_triplets[inst_uriref_relation][property][object] == 1 {
                                        	                        	set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                	                	num_ln++
									}
									biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"], biolinkobject)
                                                                        map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                                                                        if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
										set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                                                num_ln++
                                                                        }
								}
							}
						} else {
							fmt.Println("Warning: unmapped biosample ", tertiarykey)
						}
					} else if secondarykey == "object" {
						object = rdf.CompU("http://uniprot.org/uniprot/", tertiarykey)
						triplet := rdf.FormT(uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(uriref_relation, seckey2uri[secondarykey], object)
						if map_triplets[uriref_relation][seckey2uri[secondarykey]][object] == 1 {
			                        	output_file.WriteString(triplet)
							num_ln++
						}

						triplet = rdf.FormT(crm, uris_graph["tlp2tlp"], object)
						map_triplets.Add(crm, uris_graph["tlp2tlp"], object)
						output_file.WriteString(triplet)
                                                num_ln++

                                                triplet = rdf.FormT(object, uris_graph["tlp2tlp"], crm)
                                                map_triplets.Add(object, uris_graph["tlp2tlp"], crm)
                                                output_file.WriteString(triplet)
                                                num_ln++

						triplet = rdf.FormT(object, uris_graph["sub2cls"], uris_graph["tfactor"])
                                                map_triplets.Add(object, uris_graph["sub2cls"], uris_graph["tfactor"])
						if map_triplets[object][uris_graph["sub2cls"]][uris_graph["tfactor"]] == 1 {
                                                	output_file.WriteString(triplet)
                                                        num_ln++
						}

						old_object := object
						object = rdf.CompU("http://purl.uniprot.org/uniprot/", tertiarykey)
						triplet = rdf.FormT(old_object, uris_graph["sth2exm"], object)
                                                map_triplets.Add(old_object, uris_graph["sth2exm"], object)
                                                if map_triplets[old_object][uris_graph["sth2exm"]][object] == 1 {
                                                	output_file.WriteString(triplet)
                                                        num_ln++
                                                }
						biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Protein")
                                                biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"], biolinkobject)
                                                map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                                                if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
                                                	output_file.WriteString(biolinktriplet)
                                                        num_ln++
                                                }
                                        } else if secondarykey == "subject" {
						triplet := rdf.FormT(uriref_relation, seckey2uri[secondarykey], crm)
                                                map_triplets.Add(uriref_relation, seckey2uri[secondarykey], crm)
                                                output_file.WriteString(triplet)
                                                num_ln++
					} else if secondarykey == "source" {
						object = sources2uri[tertiarykey][0]
						triplet_inst := rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
						set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++
					} else if secondarykey == "crossref" {
						object = rdf.FormU(tertiarykey)
                                                triplet_inst := rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++
					}
                                }
                        }
                }
		output_file.WriteString(set_inst)
        }

	//Write and return
        rdfFile2export.Write([]byte(output_file.String()))
        output_file.Reset()

	return map_triplets, num_ln, nil
}

//Func to export CRM-TargetGene data to an RDF file (use the map3D as intermediary) --> graph crm2gene
//Input:
        //path_source: path of the file, the table with the CRM database
        //path_rdf: path where we want to save the RDF file
        //name_source: name of the database used as source
	//crm_type: type of crm sequence ("enhancer" or "silencer")
//Output:
        //map_triplets : map3D with the CRM database data
        //num_ln : number of triplets of the RDF file
        //err: nill if the execution is correct
func ExportCRM2gene2rdf(path_source, path_rdf, source_name, crm_type string)(map_triplets util.Set3D, num_ln int, err error){
	//File creation
        rdfFile2export, err := os.Create(path_rdf)
        if err != nil {
                panic(err)
        }
        defer rdfFile2export.Close()
        var output_file strings.Builder
        keys4graph := make(util.SliceSet)

        //Annotation properties
        keys4graph["Apys"] = []string{
                "sth2dfn",
                "sth2lbl",
        }

        //Object properties
	var crm2tgene, graph_name string
	if crm_type == "enhancer" { //The enhancers regulate target genes positively
		crm2tgene = "preg2targ"
		graph_name = "crm2pgene"
	} else if crm_type == "silencer" { //The silencers regulate target genes negatively
		crm2tgene = "nreg2targ"
		graph_name = "crm2ngene"
	} else {
		crm2tgene = "reg2targ"
		graph_name = "crm2ugene"
	}
        keys4graph["Opys"] = []string{
                "ins2cls",
                "sub2cls",
                "subject",
		"predicate",
		"object",
		crm2tgene,
		"sth2evd",
		"sth2mtd",
		"sth2ori",
		"sth2src",
		"biolinkcat",
		"sth2exm",
		"sth2clm",
		"seq2sample",
        }

        //Parental classes
        keys4graph["Prns"] = []string{
		"stm",
		"db",
        }

        //Header of the RDF file with the definition of the properties and the parental classes
        num_ln = 0
        header, num_ln_header := rdf.Capita(keys4graph)
        output_file.WriteString(header)
        num_ln += num_ln_header

	//URIs for graph
        uris_graph := rdf.FmtURIs(keys4graph)

	//BGW properties to Biolink properties
	proptypes := keys4graph.Keys()
	for _,proptype := range proptypes {
		if proptype == "Apys" || proptype == "Opys" {
			properties := keys4graph[proptype]
			for _,property := range properties {
				if _, ok := propbgw2propbiolink[property]; ok {
					triplet := rdf.FormT(uris_graph[property], propbgw2propbiolink[property][0], propbgw2propbiolink[property][1])
					output_file.WriteString(triplet)
					num_ln++
				}
			}
		}
	}

	//Database
	triplet := rdf.FormT(sources2uri[source_name][0], uris_graph["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
        output_file.WriteString(triplet)

	triplet = rdf.FormT(sources2uri[source_name][0], uris_graph["sub2cls"], sources2uri[source_name][2])
        output_file.WriteString(triplet)

	triplet = rdf.FormT(sources2uri[source_name][0], uris_graph["sth2lbl"], rdf.FormL(sources2uri[source_name][1]))
	output_file.WriteString(triplet)

        triplet = rdf.FormT(sources2uri[source_name][0], uris_graph["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "InformationResource"))
	output_file.WriteString(triplet)
        num_ln = num_ln + 4

	//source to RDF
        name_fields := []string{
                        "subject",
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
                        "object",
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

        skip := []int{1,2,3,4,5,6,7,8,9,10,11,13,15,16,21,22,23,24,25,26,27,28,29}
        resource, err := SeqElem2map3D(path_source, "crm2gene", name_fields, skip)
        if err != nil {
                panic(err)
        }

        //Map secondary keys of source to URI properties
        seckey2uri := map[string]string {
                "subject":	uris_graph["subject"],
                "object":	uris_graph["object"],
                "crossref":	uris_graph["sth2ori"],
                "pubmed":	uris_graph["sth2evd"],
                "method":	uris_graph["sth2mtd"],
                "source":	uris_graph["sth2src"],
		"biosample_name":       uris_graph["seq2sample"],
        }

        //Map 3D of methods
	name_fields = []string{
		"method_label",
		"efo",
		"obi",
		"eco",
		"bao",
		"ncit",
		"mi",
	}
	skip = []int{}
        methods2ids,err := EnrichFile2map3D("../../../data/enrich/methods.tsv", name_fields, skip)
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
        biosamples, err := EnrichFile2map3D("../../../data/enrich/biosamples.tsv", name_fields, skip)
        if err != nil {
                panic(err)
        }

        //Map for gene enrich
	name_fields = []string{
		"symbol",
		"ncbigene",
                "uniprot",
	}
        symbol2enrich, err := EnrichFile2map3D("../../../data/enrich/symbol2enrich.tsv", name_fields, skip)
        if err != nil {
                panic(err)
        }

        //Map of human genes of BGW
        genes_xmap := bgw.NewXmap()
        err = genes_xmap.Unmarshal("../../../data/gene/9606.json") //build a map from a json
        util.CheckE(err)

        //Set3D of the source to triplets
	map_triplets = make(util.Set3D)
	url_graph := strings.Join([]string{"http://rdf.biogateway.eu/", graph_name, "/"}, "")
	url_crm := "http://rdf.biogateway.eu/crm/9606/"
	uriref_cls := rdf.CompU(rdf.Nss["owl"], "Class")
	graph_ids := resource.Keys()
        for _,primarykey := range graph_ids {
		set2D := resource[primarykey]
		crm := rdf.CompU(url_crm, set2D["subject"].Keys()[0])

                //The relation is a class
                uriref_relation := rdf.CompU(url_graph, primarykey)
                triplet := rdf.FormT(uriref_relation, uris_graph["ins2cls"], uriref_cls)
		map_triplets.Add(uriref_relation, uris_graph["ins2cls"], uriref_cls)
                output_file.WriteString(triplet)
                num_ln++

                //The relation is a subclass
                triplet = rdf.FormT(uriref_relation, uris_graph["sub2cls"], uris_graph["stm"])
		map_triplets.Add(uriref_relation, uris_graph["sub2cls"], uris_graph["stm"])
                output_file.WriteString(triplet)
		num_ln++

                //The relation has a predicate
                triplet = rdf.FormT(uriref_relation, uris_graph["predicate"], uris_graph[crm2tgene])
                map_triplets.Add(uriref_relation, uris_graph["predicate"], uris_graph[crm2tgene])
                output_file.WriteString(triplet)
                num_ln++

		//prefLabel of the relation
		triplet = rdf.FormT(uriref_relation, uris_graph["sth2lbl"], rdf.FormL(primarykey))
		map_triplets.Add(uriref_relation, uris_graph["sth2lbl"], rdf.FormL(primarykey))
                output_file.WriteString(triplet)
                num_ln++

		//Definition of the relation
                subject := set2D["subject"].Keys()[0]
                object := strings.Split(primarykey, "/")
		definition := fmt.Sprintf("Association between cis-regulatory module %s and gene %s", subject, object[len(object)-1])
		triplet = rdf.FormT(uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
		map_triplets.Add(uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
                output_file.WriteString(triplet)
                num_ln++

		//Intance of the relation class
		database := set2D["source"].Keys()[0]
		inst_id_relation := strings.Join([]string{primarykey, "#", strings.ToLower(database)}, "")
		inst_uriref_relation := rdf.CompU(url_graph, inst_id_relation)
                map_triplets.Add(inst_uriref_relation, uris_graph["ins2cls"], uriref_relation)
                inst2class := rdf.FormT(inst_uriref_relation, uris_graph["ins2cls"], uriref_relation)
                set_inst := fmt.Sprintf("%s", inst2class)
		num_ln++

		//prefLabel of instance
		triplet = rdf.FormT(inst_uriref_relation, uris_graph["sth2lbl"], rdf.FormL(inst_id_relation))
                map_triplets.Add(inst_uriref_relation, uris_graph["sth2lbl"], inst_id_relation)
		set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

                //Definition of the instance
                definition = fmt.Sprintf("%s according to %s", definition, database)
                triplet = rdf.FormT(inst_uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
                map_triplets.Add(inst_uriref_relation, uris_graph["sth2dfn"], rdf.FormL(definition))
		set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

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
						triplet_inst := rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++

                                                old_object := object
                                                object = rdf.CompU("http://www.ncbi.nlm.nih.gov/pubmed/", tertiarykey)
                                                triplet = rdf.FormT(old_object, uris_graph["sth2exm"], object)
                                                map_triplets.Add(old_object, uris_graph["sth2exm"], object)
                                                if map_triplets[old_object][uris_graph["sth2exm"]][object] == 1 {
							set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                                                        num_ln++
                                                }

						biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Article")
						biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"], biolinkobject)
                                                map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                                                if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
                                                	output_file.WriteString(biolinktriplet)
                                                        num_ln++
                                                }
					} else if secondarykey == "method" {
                                                if _, ok := methods2ids[tertiarykey]; ok {
							ontologies := methods2ids[tertiarykey].Keys()
                                                        for _,ontology := range ontologies {
								methods_ids := methods2ids[tertiarykey][ontology].Keys()
                                                                for _,method_id := range methods_ids {
                                                                        if method_id != "-" {
										object = rdf.CompU(methkey2uri[ontology], method_id)
										triplet_inst := rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
										map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
										if map_triplets[inst_uriref_relation][seckey2uri[secondarykey]][object] == 1 {
										//Because in the mapping process, different strings can use the same ID,
										//which would generate duplicate triplets
											set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
        	                        	               				num_ln++
										}
										if ontology == "eco" {
											biolinkobject := rdf.CompU(rdf.Nss["biolink"], "EvidenceType")
											biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"],biolinkobject)
                                                                                	map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                                                                                	if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
												set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                                                        	num_ln++
                                                                                	}
										}
                                                                        }
                                                                }
                                                        }
                                                } else {
							fmt.Println("Warning: unmapped method ", tertiarykey)
						}
                                        } else if secondarykey == "biosample_name" {
						if _, ok := biosamples[tertiarykey]; ok {
							biosamples_onto := biosamples[tertiarykey].Keys()
                                                        for _, onto := range biosamples_onto {
								biosample_id := biosamples[tertiarykey][onto].Keys()[0]
								if biosample_id != "-" {
									property := seckey2uri[secondarykey]
                                                                        switch{
                                                                        case onto == "CL":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Cell")
                                                                        case onto == "CLO":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "CellLine")
                                                                        case onto == "UBERON":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "GrossAnatomicalStructure")
									case onto == "BTO":
										switch{
										case biosamples[tertiarykey]["type"].Keys()[0] == "CL":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "CellLine")
										case biosamples[tertiarykey]["type"].Keys()[0] == "CT":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Cell")
										case biosamples[tertiarykey]["type"].Keys()[0] == "A":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "GrossAnatomicalStructure")
										}
									case onto == "type":
										continue
                                                                        }
									object = rdf.CompU("http://purl.obolibrary.org/obo/", biosample_id)
        	                                                        triplet_inst := rdf.FormT(inst_uriref_relation, property, object)
                        	                                        map_triplets.Add(inst_uriref_relation, property, object)
									if map_triplets[inst_uriref_relation][property][object] == 1 {
                                        	                        	set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                	                	num_ln++
									}
									biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"], biolinkobject)
                                                                        map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                                                                        if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
										set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                                                num_ln++
                                                                        }
								}
							}
						} else {
							fmt.Println("Warning: unmapped biosample ", tertiarykey)
						}
					} else if secondarykey == "object" {
						var gsymbgw string
						var inbgw bool
						if _, ok := genes_xmap.Lblg[tertiarykey]["bgwg"]; ok {
							gsymbgw = genes_xmap.Lblg[tertiarykey]["bgwg"].Keys()[0]
							inbgw = true
						} else if _, ok := genes_xmap.Syng[tertiarykey]["bgwg"]; ok {
							gsymbgw = genes_xmap.Syng[tertiarykey]["bgwg"].Keys()[0]
							inbgw = true
						} else if _, ok := genes_xmap.Ncbig[tertiarykey]["bgwg"]; ok {
							gsymbgw = genes_xmap.Ncbig[tertiarykey]["bgwg"].Keys()[0]
							inbgw = true
						} else {
							gsymbgw = tertiarykey
							inbgw = false
						}

						if inbgw {
							object = rdf.CompU("http://rdf.biogateway.eu/gene/", gsymbgw)
						} else {
							object = rdf.CompU("https://identifiers.org/hgnc.symbol:", gsymbgw)
						}
                                                triplet := rdf.FormT(uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(uriref_relation, seckey2uri[secondarykey], object)
                                             	output_file.WriteString(triplet)
                                                num_ln++

						triplet = rdf.FormT(crm, uris_graph[crm2tgene], object)
						map_triplets.Add(crm, uris_graph[crm2tgene], object)
                                                output_file.WriteString(triplet)
                                                num_ln++

						if _,ok := symbol2enrich[tertiarykey]; ok {
							old_object := object
							ncbigene := symbol2enrich[tertiarykey]["ncbigene"].Keys()[0]
							//Because each gene has only one ncbigene ID
							if ncbigene != "-" {
								object = rdf.CompU("http://identifiers.org/ncbigene/", ncbigene)
								triplet = rdf.FormT(old_object, uris_graph["sth2clm"], object)
								map_triplets.Add(old_object, uris_graph["sth2clm"], object)
                                        	                if map_triplets[old_object][uris_graph["sth2clm"]][object] == 1 {
                                                	               	output_file.WriteString(triplet)
                                                        	       	num_ln++
                                                        	}
								biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Gene")
        	                                          	biolinktriplet := rdf.FormT(object, uris_graph["biolinkcat"], biolinkobject)
                	                                        map_triplets.Add(object, uris_graph["biolinkcat"], biolinkobject)
                        	                                if map_triplets[object][uris_graph["biolinkcat"]][biolinkobject] == 1 {
                                	                        	output_file.WriteString(biolinktriplet)
                                        	                       	num_ln++
								}
                                                	}
						}
                                        } else if secondarykey == "subject" {
						triplet := rdf.FormT(uriref_relation, seckey2uri[secondarykey], crm)
                                                map_triplets.Add(uriref_relation, seckey2uri[secondarykey], crm)
                                                output_file.WriteString(triplet)
                                                num_ln++
					} else if secondarykey == "source" {
						object = sources2uri[tertiarykey][0]
                                                triplet = rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                                                num_ln++
					} else if secondarykey == "crossref" {
						object = rdf.FormU(tertiarykey)
                                                triplet = rdf.FormT(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                map_triplets.Add(inst_uriref_relation, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                                                num_ln++
					}
                                }
                        }
                }
		output_file.WriteString(set_inst)
        }

	//Write and return
        rdfFile2export.Write([]byte(output_file.String()))
        output_file.Reset()

	return map_triplets, num_ln, nil
}

//Func to export gene coordinates to an RDF file (use the set3D as intermediary)
//Input:
        //path_source : path of the table file with the information of the genes
        //path_rdf : path where we want to save the RDF file
        //txid : ID of the taxon. Example: "9606"
        //p_genes_xmap : pointer of bgw genes map from json file
        //p_chr_map : pointer of chromosomes map
//Output:
        //map_triplets : map3D with the triplets generated for the RDF file
        //num_ln : number of triplets of the RDF file
        //err: nill if the execution is correct
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
		"sth2clm",
		"seq2strand",
        }

        //Annotation properties
        keys4genecoord["Apys"] = []string{
                "sth2lbl",
		"sth2syn",
        }

	//Parental classes
        keys4genecoord["Prns"] = []string{
                "assembly",     //genome assembly
                "chr",          //chromosomes
                "forwstrand",	//Positive strand
                "revstrand",	//Negative strand
        }

        //URIS
        uris := rdf.FmtURIs(keys4genecoord)

	//Header of the RDF file
        num_ln = 0
        header, num_ln_header := rdf.Capita(keys4genecoord)
        output_file.WriteString(header)
        num_ln += num_ln_header

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

	//Datatype properties (the function rdf.Capita does not include this type of properties)
        datatype_prop := []string{
                "seq2start",
	        "seq2end",
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

	//Assembly
	path_assembly := "https://www.ncbi.nlm.nih.gov/assembly/"
	assembly_urirefseq := chr_map[txid][chr_keys[0]]["assembly_urirefseq"].Keys()[0]
	subject := rdf.CompU(path_assembly, assembly_urirefseq)
	triplet := rdf.FormT(subject, uris["ins2cls"], rdf.CompU(rdf.Nss["owl"], "Class"))
        output_file.WriteString(triplet)

	triplet = rdf.FormT(subject, uris["sub2cls"], "<http://purl.obolibrary.org/obo/SO_0001505>")
        output_file.WriteString(triplet)

	triplet = rdf.FormT(subject, uris["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "Genome"))
        output_file.WriteString(triplet)

	assembly_prefLabel := chr_map[txid][chr_keys[0]]["assembly_label"].Keys()[0]
        triplet = rdf.FormT(subject, uris["sth2lbl"], rdf.FormL(assembly_prefLabel))
        output_file.WriteString(triplet)

	assembly_altLabel := chr_map[txid][chr_keys[0]]["assembly_altlabel"].Keys()[0]
	if assembly_altLabel != "-" {
		triplet = rdf.FormT(subject, uris["sth2syn"], rdf.FormL(assembly_altLabel))
		output_file.WriteString(triplet)
	        num_ln = num_ln + 5
	} else {num_ln = num_ln + 4}

        //Map secondary keys of source to URI properties
        seckey2uri := map[string]string {
                "chr":          uris["seq2chr"],
                "start":        Dpys["seq2start"][0],
                "end":          Dpys["seq2end"][0],
		"assembly":	uris["seq2version"],
		"strand":	uris["seq2strand"],
        }

	//Gene coordinates to triplets
        name_fields := []string{
		"ID",
                "symbol",
		"geneid",
                "chr",
                "start",
                "end",
                "strand",
		"assembly",
        }
        skip := []int{0}
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
		var gsymbgw string
		var inbgw bool
                if _, ok := genes_xmap.Lblg[gene_symbol]["bgwg"]; ok {
                	gsymbgw = genes_xmap.Lblg[gene_symbol]["bgwg"].Keys()[0]
			inbgw = true
                } else if _, ok := genes_xmap.Syng[gene_symbol]["bgwg"]; ok {
                        gsymbgw = genes_xmap.Syng[gene_symbol]["bgwg"].Keys()[0]
			inbgw = true
                } else if _, ok := genes_xmap.Ncbig[geneid]["bgwg"]; ok {
                        gsymbgw = genes_xmap.Ncbig[geneid]["bgwg"].Keys()[0]
			inbgw = true
                } else {
			gsymbgw = gene_symbol
			inbgw = false
                       continue}

		var subject, object, predicate, triplet string
		if inbgw {
			subject = rdf.CompU("http://rdf.biogateway.eu/gene/", gsymbgw)
		} else {
			subject = rdf.CompU("https://identifiers.org/hgnc.symbol:", gsymbgw)
		}
		object = rdf.CompU("http://identifiers.org/ncbigene/", geneid)
		triplet = rdf.FormT(subject, uris["sth2clm"], object)
		map_triplets.Add(object, uris["sth2clm"], object)
		output_file.WriteString(triplet)

                biolinkobject := rdf.CompU(rdf.Nss["biolink"], "Gene")
                biolinktriplet := rdf.FormT(object, uris["biolinkcat"], biolinkobject)
                map_triplets.Add(object, uris["biolinkcat"], biolinkobject)
		output_file.WriteString(biolinktriplet)
		num_ln = num_ln + 2

		secondarykeys := set2D.Keys()
                for _,secondarykey := range secondarykeys {
			tertiarykeys := set2D[secondarykey].Keys()
                       	for _,tertiarykey := range tertiarykeys {
        	                if secondarykey == "start" || secondarykey == "end" {
                                	object = strings.Join([]string{rdf.FormL(tertiarykey),"^^<http://www.w3.org/2001/XMLSchema#integer>"},"")
                                } else if secondarykey == "chr" {
					if _, ok := chr_map[txid][tertiarykey]["chr_urirefseq"]; ok {
	                                	for key,_ := range chr_map[txid][tertiarykey]["chr_urirefseq"]{
        	                                	object = strings.Join([]string{"<https://www.ncbi.nlm.nih.gov/nuccore/", key, ">"}, "")
                	                        }
					} else {object = strings.Join([]string{"<https://www.ncbi.nlm.nih.gov/nuccore/", tertiarykey, ">"}, "")}
                                } else if secondarykey == "assembly" {
					object = rdf.CompU("https://www.ncbi.nlm.nih.gov/assembly/", tertiarykey)
				} else if secondarykey == "strand" {
					if tertiarykey == "+" {
						object = "<http://biohackathon.org/resource/faldo#ForwardStrandPosition>"
					} else if tertiarykey == "-" {
						object = "<http://biohackathon.org/resource/faldo#ReverseStrandPosition>"
					} else {continue}
				} else {continue}
                                predicate = seckey2uri[secondarykey]
                                triplet = rdf.FormT(subject, predicate, object)
                                map_triplets.Add(subject, predicate, object)
                                output_file.WriteString(triplet)
                                num_ln++
                        }
                  }
	}

	//write and return
        rdfFile2export.Write([]byte(output_file.String()))
        output_file.Reset()

        return map_triplets, num_ln, nil
}

//Func to export a TAD database to an RDF file (use the map3D as intermediary)
//Input:
        //path_source: path of the file, the table with the TAD database
        //path_rdf: path where we want to save the RDF file
        //name_source: name of the database used as source

//Output:
        //map_triplets : map3D with the TAD database data
        //num_ln : number of triplets of the RDF file
        //err: nill if the execution is correct
func ExportTAD2rdf(path_source, path_rdf, source_name string) (map_triplets util.Set3D, num_ln int, err error){
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
                "seq2txn",
                "seq2chr",
                "seq2version",
		"seq2sample",
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

                triplet = rdf.FormT(assembly2uri[value][0], uris_tad["sub2cls"], assembly2uri[value][3])
                output_file.WriteString(triplet)

                triplet = rdf.FormT(assembly2uri[value][0], uris_tad["sth2lbl"], rdf.FormL(assembly2uri[value][1]))
                output_file.WriteString(triplet)

		triplet = rdf.FormT(assembly2uri[value][0], uris_tad["sth2syn"], rdf.FormL(assembly2uri[value][2]))
                output_file.WriteString(triplet)


                triplet = rdf.FormT(assembly2uri[value][0], uris_tad["biolinkcat"], rdf.CompU(rdf.Nss["biolink"], "Genome"))
                output_file.WriteString(triplet)
                num_ln = num_ln + 5
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
        name_fields := []string{
                        "TAD_ID",
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
        skip := []int{0,1,2,3,4,9,10,11,16,17}
        resource, err := SeqElem2map3D(path_source, "tad", name_fields, skip)
        if err != nil {
                panic(err)
	}

	//Map secondary keys of source to URI properties
        seckey2uri := map[string]string {
                "chr":          uris_tad["seq2chr"],
                "start":        Dpys["seq2start"][0],
                "end":          Dpys["seq2end"][0],
                "assembly":     uris_tad["seq2version"],
                "crossref":     uris_tad["sth2ori"],
                "pubmed":       uris_tad["sth2evd"],
                "method":       uris_tad["sth2mtd"],
                "source":       uris_tad["sth2src"],
		"biosample_name":	uris_tad["seq2sample"],
        }

        //Map 3D of methods
        name_fields = []string{
                "method_label",
                "efo",
                "obi",
                "eco",
                "bao",
                "ncit",
                "mi",
        }
        skip = []int{}
        methods2ids,err := EnrichFile2map3D("../../../data/enrich/methods.tsv", name_fields, skip)
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
        biosamples, err := EnrichFile2map3D("../../../data/enrich/biosamples.tsv", name_fields, skip)
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
                chr := set2D["chr"].Keys()[0]
		start := set2D["start"].Keys()[0]
		end := set2D["end"].Keys()[0]
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
		source := set2D["source"].Keys()[0]
                instid_uriref := fmt.Sprintf("<%s%s#%s>", str_tad9606, primarykey, strings.ToLower(source))
                map_triplets.Add(instid_uriref, uris_tad["ins2cls"], tadid_uriref)
                inst2class := rdf.FormT(instid_uriref, uris_tad["ins2cls"], tadid_uriref)
                set_inst := fmt.Sprintf("%s", inst2class)
                num_ln++

                //Definition of the instances
                definition = fmt.Sprintf("%s according to %s", definition, source)
                triplet = rdf.FormT(instid_uriref, uris_tad["sth2dfn"], rdf.FormL(definition))
                map_triplets.Add(instid_uriref, uris_tad["sth2dfn"], rdf.FormL(definition))
                set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

                //prefLabel of instances
                triplet = rdf.FormT(instid_uriref, uris_tad["sth2lbl"], rdf.FormL(strings.Join([]string{"tad#", primarykey, "#", strings.ToLower(source)}, "")))
                map_triplets.Add(instid_uriref, uris_tad["sth2lbl"], rdf.FormL(strings.Join([]string{"tad#", primarykey, "#", strings.ToLower(source)}, "")))
                set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                num_ln++

                //The sequences belong to Homo sapiens
                humantaxon_uri := rdf.CompU(rdf.Nss["obo"], "NCBITaxon_9606")
                triplet = rdf.FormT(tadid_uriref, uris_tad["seq2txn"], humantaxon_uri)
                map_triplets.Add(tadid_uriref, uris_tad["seq2txn"], humantaxon_uri)
                output_file.WriteString(triplet)
		num_ln++

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
                                                triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++

                                                old_object := object
                                                object = rdf.CompU("http://www.ncbi.nlm.nih.gov/pubmed/", tertiarykey)
                                                triplet = rdf.FormT(old_object, uris_tad["sth2exm"], object)
                                                map_triplets.Add(old_object, uris_tad["sth2exm"], object)
                                                if map_triplets[old_object][uris_tad["sth2exm"]][object] == 1 {
							set_inst = fmt.Sprintf("%s%s", set_inst, triplet)
                                                        num_ln++
                                                }

                                                biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Article")
                                                biolinktriplet := rdf.FormT(object, uris_tad["biolinkcat"], biolinkobject)
                                                map_triplets.Add(object, uris_tad["biolinkcat"], biolinkobject)
                                                if map_triplets[object][uris_tad["biolinkcat"]][biolinkobject] == 1 {
							set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
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
                                                                                triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                                                map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                                                if map_triplets[instid_uriref][seckey2uri[secondarykey]][object] == 1 {
                                                                                //Because in the mapping process, different strings can use the same ID,
                                                                                //which would generate duplicate triplets
                                                                                        set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                                                        num_ln++
                                                                                }
										if ontology == "eco" {
											biolinkobject := rdf.CompU(rdf.Nss["biolink"], "EvidenceType")
                                                                                	biolinktriplet := rdf.FormT(object, uris_tad["biolinkcat"], biolinkobject)
                                                                                	map_triplets.Add(object, uris_tad["biolinkcat"], biolinkobject)
                                                                                	if map_triplets[object][uris_tad["biolinkcat"]][biolinkobject] == 1 {
												set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                                                        	num_ln++
                                                                                	}
										}
									}
                                                                }
                                                        }
                                                } else {
							fmt.Println("Warning: unmapped method ", tertiarykey)
						}
                                                continue
                                        } else if secondarykey == "biosample_name" {
                                                if _, ok := biosamples[tertiarykey]; ok {
                                                        biosamples_onto := biosamples[tertiarykey].Keys()
                                                        for _, onto := range biosamples_onto {
                                                                biosample_id := biosamples[tertiarykey][onto].Keys()[0]
                                                                if biosample_id != "-" {
									property := seckey2uri[secondarykey]
                                                                        switch{
                                                                        case onto == "CL":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Cell")
                                                                        case onto == "CLO":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "CellLine")
                                                                        case onto == "UBERON":
										biolinkobject = rdf.CompU(rdf.Nss["biolink"], "GrossAnatomicalStructure")
									case onto == "BTO":
                                                                                switch{
                                                                                case biosamples[tertiarykey]["type"].Keys()[0] == "CL":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "CellLine")
                                                                                case biosamples[tertiarykey]["type"].Keys()[0] == "CT":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "Cell")
                                                                                case biosamples[tertiarykey]["type"].Keys()[0] == "A":
											biolinkobject = rdf.CompU(rdf.Nss["biolink"], "GrossAnatomicalStructure")
                                                                                }
                                                                        case onto == "type":
                                                                                continue
                                                                        }
									object = rdf.CompU("http://purl.obolibrary.org/obo/", biosample_id)
                                                                        triplet_inst := rdf.FormT(instid_uriref, property, object)
                                                                        map_triplets.Add(instid_uriref, property, object)
                                                                        if map_triplets[instid_uriref][property][object] == 1 {
                                                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                                                num_ln++
                                                                        }
									biolinktriplet := rdf.FormT(object, uris_tad["biolinkcat"], biolinkobject)
                                                                        map_triplets.Add(object, uris_tad["biolinkcat"], biolinkobject)
                                                                        if map_triplets[object][uris_tad["biolinkcat"]][biolinkobject] == 1 {
										set_inst = fmt.Sprintf("%s%s", set_inst, biolinktriplet)
                                                                                num_ln++
                                                                        }
                                                                }
                                                        }
                                                } else {
							fmt.Println("Warning: unmapped biosample ", tertiarykey)
						}
                                                continue
					} else if secondarykey == "assembly" {
                                                object = assembly2uri[tertiarykey][0]
                                        } else if secondarykey == "chr" {
                                                object = chromosomes[tertiarykey][0]
                                        } else if secondarykey == "source" {
                                                object = sources2uri[tertiarykey][0]
                                                triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++
                                                continue
                                        } else if secondarykey == "start" || secondarykey == "end" {
                                                object = strings.Join([]string{rdf.FormL(tertiarykey),"^^<http://www.w3.org/2001/XMLSchema#integer>"},"")
					} else if secondarykey == "crossref" {
                                        	object = rdf.FormU(tertiarykey)
                                                triplet_inst := rdf.FormT(instid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(instid_uriref, seckey2uri[secondarykey], object)
                                                set_inst = fmt.Sprintf("%s%s", set_inst, triplet_inst)
                                                num_ln++
                                                continue
					} else {
                                                object = rdf.FormL(tertiarykey)
                                        }

					if object != "-" {
                                                triplet := rdf.FormT(tadid_uriref, seckey2uri[secondarykey], object)
                                                map_triplets.Add(tadid_uriref, seckey2uri[secondarykey], object)
                                                output_file.WriteString(triplet)
						num_ln++

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

        return map_triplets, num_ln, nil
}

func main(){
	//Set of CRM sources
	enh_sources := util.SliceSet{
		"CancerEnD":	{"../../../databases/enh/CancerEnD.tsv", "../../../rdfs/enh/CancerEnD.nt"},
                "ChromHMM":     {"../../../databases/enh/ChromHMM.tsv", "../../../rdfs/enh/ChromHMM.nt"},
		"DiseaseEnhancer":	{"../../../databases/enh/DiseaseEnhancer.tsv", "../../../rdfs/enh/DiseaseEnhancer.nt"},
		"ENdb":		{"../../../databases/enh/ENdb.tsv", "../../../rdfs/enh/ENdb.nt"},
		"EnDisease":	{"../../../databases/enh/EnDisease.tsv", "../../../rdfs/enh/EnDisease.nt"},
		"EnhFFL":	{"../../../databases/enh/EnhFFL.tsv", "../../../rdfs/enh/EnhFFL.nt"},
                "EnhancerAtlas":        {"../../../databases/enh/EnhancerAtlas.tsv", "../../../rdfs/enh/EnhancerAtlas.nt"},
		"EnhancerDB":	{"../../../databases/enh/EnhancerDB.tsv", "../../../rdfs/enh/EnhancerDB.nt"},
		"Ensembl":	{"../../../databases/enh/Ensembl.tsv", "../../../rdfs/enh/Ensembl.nt"},
		"FANTOM5":	{"../../../databases/enh/FANTOM5.tsv", "../../../rdfs/enh/FANTOM5.nt"},
                "FOCS": {"../../../databases/enh/FOCS.tsv", "../../../rdfs/enh/FOCS.nt"},
		"GeneHancer":   {"../../../databases/enh/GeneHancer.tsv", "../../../rdfs/enh/GeneHancer.nt"},
                "GenoSTAN":     {"../../../databases/enh/GenoSTAN.tsv", "../../../rdfs/enh/GenoSTAN.nt"},
                "HACER":        {"../../../databases/enh/HACER.tsv", "../../../rdfs/enh/HACER.nt"},
		"JEME": {"../../../databases/enh/JEME.tsv", "../../../rdfs/enh/JEME.nt"},
		"RAEdb":	{"../../../databases/enh/RAEdb.tsv", "../../../rdfs/enh/RAEdb.nt"},
		"RefSeq":	{"../../../databases/enh/RefSeq.tsv", "../../../rdfs/enh/RefSeq.nt"},
		"Roadmap":	{"../../../databases/enh/Roadmap.tsv", "../../../rdfs/enh/Roadmap.nt"},
		"SCREEN":	{"../../../databases/enh/SCREEN.tsv", "../../../rdfs/enh/SCREEN.nt"},
                "SEA":  {"../../../databases/enh/SEA.tsv", "../../../rdfs/enh/SEA.nt"},
		"SEdb":	{"../../../databases/enh/SEdb.tsv", "../../../rdfs/enh/SEdb.nt"},
                "TiED": {"../../../databases/enh/TiED.tsv", "../../../rdfs/enh/TiED.nt"},
		"VISTA":	{"../../../databases/enh/VISTA.tsv", "../../../rdfs/enh/VISTA.nt"},
		"dbSUPER":	{"../../../databases/enh/dbSUPER.tsv", "../../../rdfs/enh/dbSUPER.nt"},
		"scEnhancer":	{"../../../databases/enh/scEnhancer.tsv", "../../../rdfs/enh/scEnhancer.nt"},
	}

	//CRM tables to RDF
        keys_enh_sources := enh_sources.Keys()
        for _, enh_source := range keys_enh_sources {
                paths := enh_sources[enh_source]
                _, num_ln, err := ExportCRM2rdf(paths[0], paths[1], enh_source)
                if err != nil {
                        panic(err)
                }
                fmt.Printf("Number of enh/CRM triplets in %v: %v\n", enh_source, num_ln)
        }

	//Relations sources
	crm2phen_sources := util.SliceSet{
        	"ENdb":         {"../../../databases/enh/ENdb.tsv", "../../../rdfs/crm2phen/ENdb.nt"},
        	"EnDisease":    {"../../../databases/enh/EnDisease.tsv", "../../../rdfs/crm2phen/EnDisease.nt"},
		"DiseaseEnhancer":	{"../../../databases/enh/DiseaseEnhancer.tsv", "../../../rdfs/crm2phen/DiseaseEnhancer.nt"},
	}

	keys_crm2phen_sources := crm2phen_sources.Keys()
        for _, crm2phen_source := range keys_crm2phen_sources {
                paths := crm2phen_sources[crm2phen_source]
                _, num_ln, err := ExportCRM2phen2rdf(paths[0], paths[1], crm2phen_source)
                if err != nil {
                        panic(err)
                }
                fmt.Printf("Number of crm2phen triplets in %v: %v\n", crm2phen_source, num_ln)
        }

	crm2tfac_sources := util.SliceSet{
		"ENdb":         {"../../../databases/enh/ENdb.tsv", "../../../rdfs/crm2tfac/ENdb.nt"},
		"EnhFFL":       {"../../../databases/enh/EnhFFL.tsv", "../../../rdfs/crm2tfac/EnhFFL.nt"},
	}

	keys_crm2tfac_sources := crm2tfac_sources.Keys()
        for _, crm2tfac_source := range keys_crm2tfac_sources {
                paths := crm2tfac_sources[crm2tfac_source]
                _, num_ln, err := ExportCRM2tfac2rdf(paths[0], paths[1], crm2tfac_source)
                if err != nil {
                        panic(err)
                }
                fmt.Printf("Number of crm2tfac triplets in %v: %v\n", crm2tfac_source, num_ln)
        }

	crm2pgene_sources := util.SliceSet{
		"CancerEnD":	{"../../../databases/enh/CancerEnD.tsv", "../../../rdfs/crm2pgene/CancerEnD.nt"},
		"DiseaseEnhancer":	{"../../../databases/enh/DiseaseEnhancer.tsv", "../../../rdfs/crm2pgene/DiseaseEnhancer.nt"},
		"ENdb":		{"../../../databases/enh/ENdb.tsv", "../../../rdfs/crm2pgene/ENdb.nt"},
		"EnDisease":	{"../../../databases/enh/EnDisease.tsv", "../../../rdfs/crm2pgene/EnDisease.nt"},
		"EnhFFL":	{"../../../databases/enh/EnhFFL.tsv", "../../../rdfs/crm2pgene/EnhFFL.nt"},
                "EnhancerAtlas":        {"../../../databases/enh/EnhancerAtlas.tsv", "../../../rdfs/crm2pgene/EnhancerAtlas.nt"},
		"FANTOM5":	{"../../../databases/enh/FANTOM5.tsv", "../../../rdfs/crm2pgene/FANTOM5.nt"},
                "FOCS": {"../../../databases/enh/FOCS.tsv", "../../../rdfs/crm2pgene/FOCS.nt"},
		"GeneHancer":   {"../../../databases/enh/GeneHancer.tsv", "../../../rdfs/crm2pgene/GeneHancer.nt"},
                "HACER":        {"../../../databases/enh/HACER.tsv", "../../../rdfs/crm2pgene/HACER.nt"},
		"JEME":	{"../../../databases/enh/JEME.tsv", "../../../rdfs/crm2pgene/JEME.nt"},
                "SEA":  {"../../../databases/enh/SEA.tsv", "../../../rdfs/crm2pgene/SEA.nt"},
		"SEdb":	{"../../../databases/enh/SEdb.tsv", "../../../rdfs/crm2pgene/SEdb.nt"},
		"VISTA":	{"../../../databases/enh/VISTA.tsv", "../../../rdfs/crm2pgene/VISTA.nt"},
		"dbSUPER":	{"../../../databases/enh/dbSUPER.tsv", "../../../rdfs/crm2pgene/dbSUPER.nt"},
		"scEnhancer":	{"../../../databases/enh/scEnhancer.tsv", "../../../rdfs/crm2pgene/scEnhancer.nt"},
	}

	keys_crm2pgene_sources := crm2pgene_sources.Keys()
        for _, crm2pgene_source := range keys_crm2pgene_sources {
                paths := crm2pgene_sources[crm2pgene_source]
                _, num_ln, err := ExportCRM2gene2rdf(paths[0], paths[1], crm2pgene_source, "enhancer")
                if err != nil {
                        panic(err)
                }
                fmt.Printf("Number of crm2pgene triplets in %v: %v\n", crm2pgene_source, num_ln)
        }

	//Set of TAD sources
	tad_sources := util.SliceSet{
		"3DGB":	{"../../../databases/tad/3DGB.tsv", "../../../rdfs/tad/3DGB.nt"},
		"TADKB":	{"../../../databases/tad/TADKB.tsv", "../../../rdfs/tad/TADKB.nt"},
	}

	//TAD tables to RDF
        keys_tad_sources := tad_sources.Keys()
        for _, tad_source := range keys_tad_sources {
                paths := tad_sources[tad_source]
		_, num_ln, err := ExportTAD2rdf(paths[0], paths[1], tad_source)
		if err != nil {
                        panic(err)
                }
                fmt.Printf("Number of TAD triplets in %v: %v\n", tad_source, num_ln)
        }

	//Generate chr map
	map_chr, err := Chr2map("../../../data/gene/chromosomes.tsv")
	if err != nil {
        	panic(err)
        }

	//Genes coordinates
	files, err := ioutil.ReadDir("../../../databases/gene/")
	if err != nil {
        	panic(err)
	}
    	for _, file := range files {
		txid := strings.Split(file.Name(), "_")[1]
		txid = strings.Split(txid, ".")[0]

		genes_xmap := bgw.NewXmap()
        	err := genes_xmap.Unmarshal(strings.Join([]string{"../../../mappings/", txid, ".json"}, ""))
        	util.CheckE(err)

		_, num_ln, err := ExportGeneCoord2rdf(strings.Join([]string{"../../../databases/gene/", file.Name()}, ""),
					strings.Join([]string{"../../../rdfs/gene/gene_coord_", txid, ".nt"}, ""),
					txid, &genes_xmap, &map_chr)

		if err != nil {
                	panic(err)
        	}
      		fmt.Printf("Number of triplets in Gene coordinates of txid %v: %v\n", txid, num_ln)
	}
}
