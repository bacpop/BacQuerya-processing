# BacQuerya

BacQuerya is a search engine that aims to consolidate and present all publicly available genomic metadata for bacterial pathogens. BacQuerya is currently in beta and as such, is unstable in some circumstances and only houses *S. pneumoniae* genomic metadata at this time. 

## BacQuerya-processing

The BacQuerya-processing pipeline sources genomic metadata by accession ID from public repositories that include: NCBI GenBank, BioSample, the European Nucleotide Archive and the Sequence Read Archive. Metadata from each of these sources is extracted and combined into a JSON document, that is then indexed using Elastic cloud ([https://www.elastic.co](https://www.elastic.co)) and is searchable from the BacQuerya website ([https://github.com/bacpop/BacQuerya](https://github.com/bacpop/BacQuerya)). 

