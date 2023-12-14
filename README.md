# BioGateway for gene regulation knowledge

## Introduction
BioGateway is a knowledge network based on graphs and RDF triples to query biological information corresponding to *Homo sapiens* and other organisms. BioGateway is a model that uses semantic web technologies to integrate information from different databases through a model that, by reusing available resources, aims to standardize biological data and improve web interoperability. In addition, the knowledge network includes biological coordinates to allow the development of query strategies that exploit the location of sequences, an aspect that other semantic biological knowledge networks do not usually exploit.

The knowledge network is available through its endpoint [https://biogateway.eu/sparql-endpoint/](https://biogateway.eu/sparql-endpoint/) and SPARQL. A short SPARQL [tutorial](./SPARQL_Tutorial.pdf) is available to introduce potencial users to this query language as well as to the use of the knowledge network.

With this contribution we also encourage the community to work on the development of different interoperable knowledge networks to connect domains and allow the development of complex federated queries.

## Structure and content

The generated knowledge network is structured in graphs, being each graph a different information domain. We distinguish two types of graphs in the network: entity graphs and relation graphs. The first ones aim to model different biological entities, while the second ones model relations between different entities.

![Graph_types](./images/graphs.png)

The knowledge network has the following [graphs](./hBGW_graphs.xlsx):
- crm : Cis Regulatory Modules. Currently only enhancer sequences, that increase gene transcription levels.
- crm2phen: Relations between CRM and phenotypes.
- crm2gene: Relations between CRM and target genes.
- crm2tfac: Relations between CRM and transcription factors.
- tad : Topologically associating domain. Domains of genome structure and regulation.
- gene : Genes.
- prot : Proteins.
- omim : OMIM ontology (phenotypes, among others).
- go : GO ontology (biological processes, molecular functions and cellular components).
- mi : Molecular Interaction Ontology.
- taxon : NCBI Taxon Ontology.
- gene2phe : Genes - Phenotypes (omim) relations .
- tfac2gene : Relations between Transcription factors (TF) and their target genes.
- prot2prot : Protein-protein interactions.
- reg2targ : Protein - Protein regulatory relations
- prot2cc : Protein - Celullar components relations.
- prot2bp : Protein - Biological processes relations.
- prot2mf : Protein - Molecular functions relations.
- ortho : Protein-protein orthology relations.  

Network statistics can also be found [here](./hBGW_statistics.xlsx).

BioGateway is an integrative model in which classes represent entities and instances represent individuals of the classes. Classes have relations with those entities and data that model the entity. Instances collect the information from each original database, so they have the relations specific of each individual. Therefore, they are useful nodes for representing the information from each source and its metadata.

## Workflow

The current version of BioGateway integrates domains associated with gene regulation (CRM and TADs). It also extends the gene graph to include aspects of genome location and provides relations to the [Biolink model](https://biolink.github.io/biolink-model/), focused on the standardization of knowledge graphs.

The following steps were taken to generate the new hBGW domains:
1. Compilation of the [databases used as sources](./images/databases.PNG).
2. Preprocessing of files [(Rcodes folder)](./Rcodes).
3. Generation of RDF files [(cisreg folder)](./cisreg). This code requires functions of the project: https://github.com/vlmir/bgw3. An RDF file fragment corresponding to [enhancers](./RDF_file_examples/CRM.nt) and [TADs](./RDF_file_examples/TAD.nt) is included in the repository as an illustrative example of output.
4. Loading and generation of graphs in [Virtuoso](https://github.com/openlink/virtuoso-opensource).

![workflow](./images/workflow.PNG)

## Use Cases

We demonstrate the use of the knowledge network in a set of advanced queries (three use cases) that were not possible in an efficient way until now. Below we provide the files corresponding to the queries and the results obtained:

- Use Case 1: [Queries](./Use_Cases/UC1/Queries.txt) - [Results](./Use_Cases/UC1/Results.xlsx).
- Use Case 2: [Queries](./Use_Cases/UC2/Queries.txt) - [Results](./Use_Cases/UC2/Results.xlsx).
- Use Case 3: [Queries](./Use_Cases/UC3/Queries.txt) - [Results](./Use_Cases/UC3/Results.tsv).

![UseCases](./images/UseCases.PNG)
