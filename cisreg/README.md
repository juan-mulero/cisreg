Running the go code requires the installation of [go](https://go.dev/) language and functions from the https://github.com/vlmir/bgw3.git repository. 

To clone the repository and execute the code:

```
git clone https://github.com/juan-mulero/cisreg.git
cd cisreg/cisreg/
go run cisreg.go
```
Execution requires to specify input and output paths for files. Also the paths of the files used for the enrichment and mapping of data. Therefore, check these paths before execution.

Recommendation: to avoid memory saturation, split databases by chromosomes if the databases to be converted to RDF are large.
