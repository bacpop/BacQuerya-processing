# BacQuerya

BacQuerya is a search engine that aims to consolidate and present all publicly available genomic metadata for bacterial pathogens. BacQuerya is currently in beta and as such, is unstable in some circumstances and only houses *S. pneumoniae* genomic metadata at this time. 

## BacQuerya-processing

The BacQuerya-processing pipeline sources genomic metadata by accession ID from public repositories that include: NCBI GenBank, BioSample, the European Nucleotide Archive and the Sequence Read Archive. Metadata from each of these sources is extracted and combined into a JSON document, that is then indexed using Elastic cloud ([https://www.elastic.co](https://www.elastic.co)) and is searchable from the BacQuerya website ([https://github.com/bacpop/BacQuerya](https://github.com/bacpop/BacQuerya)). 

Scripts in this repository may be used either with the included snakemake pipeline or individually using the helper scripts.

## Installation

To download BacQuerya-processing from source, run:
```git clone https://github.com/bacpop/BacQuerya-processing.git```

## Snakemake pipeline

Parameters for the automated Snakemake pipeline can be adjusted by modifying the ```config.yml``` file, or from the command line. An example command to run the ```retrieve_ena_read_metadata``` rule on 7 threads would be:

```snakemake --cores 1 --config n_cpu=7 retrieve_ena_read_metadata```

## Scripts

Snakemake rule:
```
retrieve_assembly_stats
```
Equivalent script command:
```
python extract_entrez_information-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}
```

Snakemake rule:
```
retrieve_annotations
```
Equivalent script command:
```
python extract_entrez_information-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}
```

Snakemake rule:
```
retrieve_genomes
```
Equivalent script command:
```
python extract_entrez_information-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}
```

