Running the go code requires the installation of [go](https://go.dev/) language and functions from the https://github.com/vlmir/bgw3.git repository. 

To clone the repository and execute the code, include cisreg as a submodule of bgw3:

```
git clone https://github.com/vlmir/bgw3.git
git clone https://github.com/juan-mulero/cisreg.git
mv cisreg/cisreg bgw3/src/
cd bgw3/src/cisreg
go run cisreg.go
```
Execution requires to specify input and output paths for files. Also the paths of the files used for the enrichment and mapping of data. Therefore, check these paths before execution.

The tdata folder contains example files, while cisreg_test.go allows to test the running of the code provided.