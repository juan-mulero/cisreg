package main

import (
        "fmt"
        "testing"
        "github.com/vlmir/bgw3/src/bgw"
        "github.com/vlmir/bgw3/src/util"
)

func Test_Enh2map3D(t *testing.T){
	//Set example types
        type set_example struct{
                primkey, seckey_chr, seckey_disease, seckey_target_gene string
                tertkey_chr, tertkey_disease, tertkey_target_gene string
                value_chr, value_disease, value_target_gene int
        }
	//Set example values
        set_ex := set_example{"CRMHS00000000001", "chr", "disease", "target_gene",
                                "chr1", "DOID:1936", "VCAM1",
                                1, 1, 1}
	//Map old2new_IDs
        name_fields := []string{
                "old_ID",
                "new_ID",
        }
        var skip []int
        old2new_IDs, err := EnrichFile2map3D("/home/jmulero/bgw3/src/cisreg/data/crm/crm_old2new_test.tsv", name_fields, skip)

	//Map new2old_IDs
        name_fields = []string{
                "new_ID",
                "old_ID",
        }
        new2old_IDs, err := EnrichFile2map3D("/home/jmulero/bgw3/src/cisreg/data/crm/crm_new2old_test.tsv", name_fields, skip)


	//Enhancers to map3D
        map3D, _, _, err := Enh2map3D("/home/jmulero/bgw3/src/cisreg/data/crm/enh_testfile.tsv", &old2new_IDs, &new2old_IDs)

	//Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }

        if map3D[set_ex.primkey][set_ex.seckey_chr][set_ex.tertkey_chr] != set_ex.value_chr{
                t.Error("Unexpected value in the map with primkey:", set_ex.primkey, ", seckey:", set_ex.seckey_chr, "and tertkey:", set_ex.tertkey_chr)
        }

        if map3D[set_ex.primkey][set_ex.seckey_disease][set_ex.tertkey_disease] != set_ex.value_disease{
                t.Error("Unexpected value in the map with primkey:", set_ex.primkey, ", seckey:", set_ex.seckey_disease,
                        "and tertkey:", set_ex.tertkey_disease)
        }
        if map3D[set_ex.primkey][set_ex.seckey_target_gene][set_ex.tertkey_target_gene] != set_ex.value_target_gene{
                t.Error("Unexpected value in the map with primkey:", set_ex.primkey, ", seckey:", set_ex.seckey_target_gene,
                        "and tertkey:", set_ex.tertkey_target_gene)
        }
}

func Test_TAD2map3D(t *testing.T){
	//Set example types
        type set_example struct{
                primarykey, secondarykey, tertiarykey string
                value int
        }

	//Set example values
        set_ex := set_example{"TADHS00000000001", "pubmed", "30871473", 1}

        //Map old2new_IDs
        name_fields := []string{
                "old_ID",
                "new_ID",
        }
        var skip []int
        old2new_IDs, err := EnrichFile2map3D("/home/jmulero/bgw3/src/cisreg/data/tad/tad_old2new_test.tsv", name_fields, skip)

        //Map new2old_IDs
        name_fields = []string{
                "new_ID",
                "old_ID",
        }
        new2old_IDs, err := EnrichFile2map3D("/home/jmulero/bgw3/src/cisreg/data/tad/tad_new2old_test.tsv", name_fields, skip)

	//TAD to map3D
        map3D, _, _, err := TAD2map3D("/home/jmulero/bgw3/src/cisreg/data/tad/tad_testfile.tsv", &old2new_IDs, &new2old_IDs)

	//Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }
        if map3D[set_ex.primarykey][set_ex.secondarykey][set_ex.tertiarykey] != set_ex.value{
                t.Error("Unexpected value in the map with primkey:", set_ex.primarykey, ", seckey:", set_ex.secondarykey,
                        "and tertkey:", set_ex.tertiarykey)
        }
}

func Test_EnrichFile2map3D(t *testing.T){
	//Set example types
        type set_example struct{
                primarykey, secondarykey, tertiarykey string
                value int
        }

	//Set example values
        set_ex := set_example{"4C", "obi", "OBI_0002458", 1}

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
        map3D,err := EnrichFile2map3D("/home/jmulero/bgw3/src/cisreg/data/enrich/methods_test.tsv", name_fields, skip)

	//Map test
        if err != nil{
                t.Error(fmt.Sprintf("Error, the map is wrong: %v", err))
        }
        if map3D[set_ex.primarykey][set_ex.secondarykey][set_ex.tertiarykey] != set_ex.value{
                t.Error("Unexpected value in the map with primarykey:", set_ex.primarykey, ", secondarykey:", set_ex.secondarykey,
                        "and tertiarykey:", set_ex.tertiarykey)
        }
}

func Test_ExportEnh2rdf(t *testing.T){
	//Set example types
	type set_example struct{
		subject, predicate, object string
		value int
	}

        //Set example values
	set_ex := set_example{
                "<http://rdf.biogateway.eu/crm/9606/CRMHS00000000001>",
                "<http://purl.obolibrary.org/obo/RO_0000052>",
                "<http://purl.obolibrary.org/obo/NCBITaxon_9606>",
                1,
        }

	//Map
	genes_xmap := bgw.NewXmap()
        err := genes_xmap.Unmarshal("/home/jmulero/bgw3/src/cisreg/data/gene/9606_test.json")
        util.CheckE(err)

        name_fields := []string{
                "old_ID",
                "new_ID",
        }
        var skip []int
        old2new_IDs, err := EnrichFile2map3D("/home/jmulero/bgw3/src/cisreg/data/crm/crm_old2new_test.tsv", name_fields, skip)

        //Map new2old_IDs
        name_fields = []string{
                "new_ID",
                "old_ID",
        }
        new2old_IDs, err := EnrichFile2map3D("/home/jmulero/bgw3/src/cisreg/data/crm/crm_new2old_test.tsv", name_fields, skip)

	map_triples, _, _, _, err := ExportEnh2rdf("/home/jmulero/bgw3/src/cisreg/data/crm/enh_testfile.tsv", "/home/jmulero/bgw3/src/cisreg/data/crm/enh2rdf.nt", "ENdb", &genes_xmap, &old2new_IDs, &new2old_IDs)
	if err != nil{
                t.Error(err)
        }

	//Map test
        if map_triples[set_ex.subject][set_ex.predicate][set_ex.object] != set_ex.value{
                t.Error("Unexpected value in the map with subject:", set_ex.subject, ", predicate:", set_ex.predicate,
                        "and object:", set_ex.object)
        }
}

func Test_ExportTAD2rdf(t *testing.T){
        //Set example types
	type set_example struct{
                subject, predicate, object string
                value int
        }

	//Set example values
	set_ex := set_example{
                "<http://rdf.biogateway.eu/tad/9606/TADHS00000000001>",
		"<http://semanticscience.org/resource/SIO_000772>",
		"<http://identifiers.org/pubmed/30871473>",
                1,
        }

	//Map
	//Map old2new_IDs
        name_fields := []string{
                "old_ID",
                "new_ID",
        }
        var skip []int
        old2new_IDs, err := EnrichFile2map3D("/home/jmulero/bgw3/src/cisreg/data/tad/tad_old2new_test.tsv", name_fields, skip)

        //Map new2old_IDs
        name_fields = []string{
                "new_ID",
                "old_ID",
        }
        new2old_IDs, err := EnrichFile2map3D("/home/jmulero/bgw3/src/cisreg/data/tad/tad_new2old_test.tsv", name_fields, skip)

	map_triples, _, _, _, err := ExportTAD2rdf("/home/jmulero/bgw3/src/cisreg/data/tad/tad_testfile.tsv", "/home/jmulero/bgw3/src/cisreg/data/tad/tad2rdf.nt", "TADKB", &old2new_IDs, &new2old_IDs)

	//Map test
	if err != nil{
                t.Error(err)
        }
        if map_triples[set_ex.subject][set_ex.predicate][set_ex.object] != set_ex.value{
                t.Error("Unexpected value in the map with subject:", set_ex.subject, ", predicate:", set_ex.predicate,
                        "and object:", set_ex.object)
        }
}

func Test_ExportGeneCoord2rdf(t*testing.T){
	//Set example types
	type set_example struct{
		subject, predicate, object string
		value int
	}

	//Set example values
	set_ex := set_example{
		"<http://rdf.biogateway.eu/gene/9606/U2AF1>",
		"<http://purl.obolibrary.org/obo/OGI_1000003>",
		"43108291",
		1,
	}

	//Maps
        genes_xmap := bgw.NewXmap()
        err := genes_xmap.Unmarshal("/home/jmulero/bgw3/src/cisreg/data/gene/9606_test.json")
        util.CheckE(err)

	chr_map, err := Chr2map("/home/jmulero/bgw3/src/cisreg/data/gene/chr_test.tsv")
        if err != nil {
                panic(err)
        }

	map_triples, _, err := ExportGeneCoord2rdf("/home/jmulero/bgw3/src/cisreg/data/gene/NCBITaxon_9606.tsv",
							"/home/jmulero/bgw3/src/cisreg/data/gene/genecoord_test.nt", "9606", &genes_xmap, &chr_map)
        if err != nil{
                t.Error(err)
        }

        if map_triples[set_ex.subject][set_ex.predicate][set_ex.object] != set_ex.value{
                t.Error("Unexpected value in the map with subject:", set_ex.subject, ", predicate:", set_ex.predicate,
                        "and object:", set_ex.object)
        }
}
