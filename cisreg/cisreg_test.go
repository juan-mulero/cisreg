package main

import (
        "fmt"
        "testing"
	"github.com/vlmir/bgw3/src/util"
	"github.com/vlmir/bgw3/src/bgw"
)

func ExistInStringSlice(set []string, search string) bool{
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

func Test_SeqElem2map3D(t *testing.T){
	//Set example types
	type set_example struct{
                primkey, seckey_chr, seckey_disease, seckey_tgene string
                tertkey_chr, tertkey_disease, tertkey_tgene string
	}
	//Set example values
        set_ex := set_example{"CRMHS00000005327", "chr", "disease", "target_gene",
                                "chr1", "DOID:1936", "VCAM1"}

	//Input data
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
        skip := []int{}
        map3D, err := SeqElem2map3D("./tdata/input/crm_db.tsv", "crm", name_fields, skip)

	//Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }

	var check bool
	check = ExistInStringSlice(map3D[set_ex.primkey][set_ex.seckey_chr].Keys(), set_ex.tertkey_chr)
	if !check {
		t.Error("Unexpected value in the map with primkey:", set_ex.primkey, ", seckey:", set_ex.seckey_chr, "and tertkey:", set_ex.tertkey_chr)
	}

	check = ExistInStringSlice(map3D[set_ex.primkey][set_ex.seckey_disease].Keys(), set_ex.tertkey_disease)
	if !check {
                t.Error("Unexpected value in the map with primkey:", set_ex.primkey, ", seckey:", set_ex.seckey_disease, "and tertkey:", set_ex.tertkey_disease)
        }

	check = ExistInStringSlice(map3D[set_ex.primkey][set_ex.seckey_tgene].Keys(), set_ex.tertkey_tgene)
        if !check {
                t.Error("Unexpected value in the map with primkey:", set_ex.primkey, ", seckey:", set_ex.seckey_tgene, "and tertkey:", set_ex.tertkey_tgene)
        }
}

func Test_Chr2map(t *testing.T){
	//Set example types
	type set_example struct{
		primkey, seckey, tertkey, quatkey string
	}
	//Set example values
        set_ex := set_example{"9606", "chr6", "chr_urirefseq", "NC_000006.12"}

	//Input data
	map_chr, err := Chr2map("./tdata/input/chromosomes.tsv")

	//Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }

        var check bool
        check = ExistInStringSlice(map_chr[set_ex.primkey][set_ex.seckey][set_ex.tertkey].Keys(), set_ex.quatkey)
        if !check {
                t.Error("Unexpected value in the map with primkey:", set_ex.primkey, ", seckey:", set_ex.seckey, ", tertkey:", set_ex.tertkey, "and quatkey:", set_ex.quatkey)
        }
}

func Test_EnrichFile2map3D(t *testing.T){
        //Set example types
        type set_example struct{
                primkey, seckey, tertkey string
        }

        //Set example values
        set_ex := set_example{"4C", "obi", "OBI_0002458"}

        //Map
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
        map3D,err := EnrichFile2map3D("./tdata/input/methods.tsv", name_fields, skip)

        //Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
	}
        var check bool
        check = ExistInStringSlice(map3D[set_ex.primkey][set_ex.seckey].Keys(), set_ex.tertkey)
        if !check {
                t.Error("Unexpected value in the map with primkey:", set_ex.primkey, ", seckey:", set_ex.seckey, "and tertkey:", set_ex.tertkey)
        }
}

func Test_ExportCRM2rdf(t *testing.T){
        //Set example types
        type set_example struct{
                subject, predicate, object string
        }

        //Set example values
        set_ex := set_example{
                "<http://rdf.biogateway.eu/crm/9606/CRMHS00000005327>",
                "<http://purl.obolibrary.org/obo/RO_0002162>",
                "<http://purl.obolibrary.org/obo/NCBITaxon_9606>",
        }

	//Map
        map3D, _, err := ExportCRM2rdf("./tdata/input/crm_db.tsv", "./tdata/output/crm.nt", "ENdb")

        //Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }
        var check bool
        check = ExistInStringSlice(map3D[set_ex.subject][set_ex.predicate].Keys(), set_ex.object)
        if !check {
                t.Error("Unexpected triplet in the map with subject:",set_ex.subject,", predicate:",set_ex.predicate,"and object:",set_ex.object)
        }
}

func Test_ExportCRM2phen2rdf(t *testing.T){
	//Set example types
        type set_example struct{
                subject, predicate, object string
        }

        //Set example values
        set_ex := set_example{
                "<http://rdf.biogateway.eu/crm2phen/bgw!CRMHS00000005327--doid!1936>",
                "<http://www.w3.org/1999/02/22-rdf-syntax-ns#subject>",
                "<http://rdf.biogateway.eu/crm/9606/CRMHS00000005327>",
        }

        //Map
        map3D, _, err := ExportCRM2phen2rdf("./tdata/input/crm_db.tsv", "./tdata/output/crm2phen.nt", "ENdb")

        //Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }
        var check bool
        check = ExistInStringSlice(map3D[set_ex.subject][set_ex.predicate].Keys(), set_ex.object)
        if !check {
                t.Error("Unexpected triplet in the map with subject:",set_ex.subject,", predicate:",set_ex.predicate,"and object:",set_ex.object)
        }
}

func Test_ExportCRM2tfac2rdf(t *testing.T){
        //Set example types
        type set_example struct{
                subject, predicate, object string
        }

        //Set example values
        set_ex := set_example{
                "<http://rdf.biogateway.eu/crm2tfac/bgw!CRMHS00000005327--uniprot!P42226>",
                "<http://www.w3.org/1999/02/22-rdf-syntax-ns#object>",
                "<http://uniprot.org/uniprot/P42226>",
        }

        //Map
        map3D, _, err := ExportCRM2tfac2rdf("./tdata/input/crm_db.tsv", "./tdata/output/crm2tfac.nt", "ENdb")

        //Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }
        var check bool
        check = ExistInStringSlice(map3D[set_ex.subject][set_ex.predicate].Keys(), set_ex.object)
        if !check {
                t.Error("Unexpected triplet in the map with subject:",set_ex.subject,", predicate:",set_ex.predicate,"and object:",set_ex.object)
        }
}

func Test_ExportCRM2gene2rdf(t *testing.T){
        //Set example types
        type set_example struct{
                subject, predicate, object string
        }

        //Set example values
        set_ex := set_example{
                "<http://rdf.biogateway.eu/crm2pgene/bgw!CRMHS00000005327--hgncsymbol!9606/VCAM1>",
                "<http://www.w3.org/1999/02/22-rdf-syntax-ns#predicate>",
                "<http://purl.obolibrary.org/obo/RO_0002429>",
        }

        //Map
        map3D, _, err := ExportCRM2gene2rdf("./tdata/input/crm_db.tsv", "./tdata/output/crm2pgene.nt", "ENdb", "enhancer")

        //Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }
        var check bool
        check = ExistInStringSlice(map3D[set_ex.subject][set_ex.predicate].Keys(), set_ex.object)
        if !check {
                t.Error("Unexpected triplet in the map with subject:",set_ex.subject,", predicate:",set_ex.predicate,"and object:",set_ex.object)
        }
}

func Test_ExportTAD2rdf(t *testing.T){
        //Set example types
        type set_example struct{
                subject, predicate, object string
        }

        //Set example values
        set_ex := set_example{
                "<http://rdf.biogateway.eu/tad/9606/TADHS00000022800>",
                "<http://purl.org/dc/terms/hasVersion>",
                "<https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26>",
        }

        //Map
        map3D, _, err := ExportTAD2rdf("./tdata/input/tad_db.tsv", "./tdata/output/tad.nt", "3DGB")

        //Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }
        var check bool
        check = ExistInStringSlice(map3D[set_ex.subject][set_ex.predicate].Keys(), set_ex.object)
        if !check {
                t.Error("Unexpected triplet in the map with subject:",set_ex.subject,", predicate:",set_ex.predicate,"and object:",set_ex.object)
        }
}

func Test_ExportGeneCoord2rdf(t *testing.T){
        //Set example types
        type set_example struct{
                subject, predicate, object string
        }

        //Set example values
        set_ex := set_example{
                "<http://rdf.biogateway.eu/gene/9606/VCAM1>",
                "<http://purl.obolibrary.org/obo/GENO_0000894>",
                "\"100719742\"^^<http://www.w3.org/2001/XMLSchema#integer>",
        }

	//Map
        map_chr, err := Chr2map("./tdata/input/chromosomes.tsv")
        if err != nil {
                panic(err)
        }

	genes_xmap := bgw.NewXmap()
        err = genes_xmap.Unmarshal("./tdata/input/9606.json")
        util.CheckE(err)

        map3D, _, err := ExportGeneCoord2rdf("./tdata/input/NCBITaxon_9606.tsv", "./tdata/output/gene_coord_9606.nt", "9606", &genes_xmap, &map_chr)

        //Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }
        var check bool
        check = ExistInStringSlice(map3D[set_ex.subject][set_ex.predicate].Keys(), set_ex.object)
        if !check {
                t.Error("Unexpected triplet in the map with subject:",set_ex.subject,", predicate:",set_ex.predicate,"and object:",set_ex.object)
        }
}
