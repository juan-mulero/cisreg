Running the go code requires the installation of [go](https://go.dev/) language and functions from the https://github.com/vlmir/bgw3.git repository. 

To clone the repository and execute the code, include cisreg as a submodule of bgw3:

```
git clone https://github.com/vlmir/bgw3.git
git clone https://github.com/juan-mulero/cisreg.git
mv cisreg/cisreg bgw3/src/
cd bgw3/src/cisreg
go run cisreg.go
```
Execution requires to specify input and output paths for files in [cisreg.go](./cisreg.go):
-  The input files are the files obtained after running the [Rcodes](../Rcodes). Also the paths of the files used for the [enrichment](./mappings/symbol2enrich.tsv) and mapping of data: [biosamples](./mappings/biosamples.tsv) and [methods](./mappings/methods.tsv).
-  The output files correspond to the paths where you want to save the RDF files.
 
Therefore, check these paths before execution.

The tdata folder contains example files for testing, while cisreg_test.go allows to test the running of the code provided.

This code generates the graphs corresponding to crm, crm2gene, crm2phen, crm2tfac, and tad. It also generates gene triplets that complement the gene graph. Although these graphs can be used independently, they were designed to be integrated into the BGW knowledge network. Consequently, the following [github repository](https://github.com/BioGateway/bgw-builder) includes the guidelines for generating the rest of BGW graphs.