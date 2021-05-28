# BacQuerya

BacQuerya is a search engine that aims to consolidate and present all publicly available genomic metadata for bacterial pathogens. BacQuerya is currently in beta and as such, is unstable in some circumstances and only houses *S. pneumoniae* genomic metadata at this time. 

# BacQuerya-processing

The BacQuerya-processing pipeline sources genomic metadata by accession ID from public repositories that include: NCBI GenBank, BioSample, the European Nucleotide Archive and the Sequence Read Archive. Metadata from each of these sources is extracted and combined into a JSON document, that is then indexed using Elastic cloud ([https://www.elastic.co](https://www.elastic.co)) and is searchable from the BacQuerya website ([https://github.com/bacpop/BacQuerya](https://github.com/bacpop/BacQuerya)). 

Scripts in this repository may be used either with the included snakemake pipeline or individually using the helper scripts.

# Installation

To install BacQuerya-processing from source, run:
```
git clone https://github.com/bacpop/BacQuerya-processing.git
conda install snakemake
conda create --file=environment.yml
conda activate snakemake
```

# Snakemake pipeline

Parameters for the automated Snakemake pipeline can be adjusted by modifying the ```config.yml``` file, or from the command line. An example command to run the ```retrieve_ena_read_metadata``` rule with 1 core on 7 threads would be:

```snakemake --cores 1 --config n_cpu=7 retrieve_ena_read_metadata```

## Rules

1. ```retrieve_assembly_stats```, ```retrieve_annotations``` and ```retrieve_genomes``` rules:

Retrieves information of interest from NCBI GenBank by accession ID. 

**Inputs**:
* ```accession_file```: Filepath of a "\n" separated list of BioSample or assembly accession IDs for asssemblies available through NCBI GenBank or Refseq. 
* ```attribute```: Retrieve assembly sequences, functional annotations or assembly statistics (specified as ```genomes```, ```annotation``` or ```assembly-stats```).
* ```email```: An email address to specify for Entrez programmatic access. 
* ```output```: Output directory name for retrieved attributes. 
* ```threads```: Number of threads for retrieval. Entrez allows up to 3 queries/second without an API key and 10 queries/second with an API key. To specify an API key, see **Adding secrets**.

**Equivalent shell command**:
```
python extract_entrez_information-runner.py -s <accession_file> -e <email> --threads <threads> -o <output> -a <attribute>
```

2. ```run_prodigal``` rule:

Predicts genes in all assemblies in an input directory and outputs predicted annotations in GFF3 format.
**Inputs**:
* ```genome_dir```: Directory containing assemblies. 
* ```output```: Output directory name for predicted annotations. 
* ```threads```: Number of threads for prediction. 

3. ```reformat_annotations``` rule:

Converts publicly available functional annotations in GFF3 format to Prokka format for direct input into Panaroo.

**Inputs**:
* ```assembly_directory```: Directory of uncompressed assembly sequences.
* ```index_file```: Filepath of JSON storing index numbers (Of the form ```{"isolateIndexNo": 0, "geneIndexNo": 0, "predictedIndexNo": 0}```).
* ```output```: Output directory name for reformatted annotations. 
* ```threads```: Number of threads for reformatting. 
* ```annotation_directory```: Directory of uncompressed functional annotations in GFF3 format. If specified, existing annotations are added to the reformatted annotations. 
* ```prodigal_directory```: Directory of functional annotation files output by prodigal. If specified, prodigal-predicted annotations are added to the reformatted annotations. If ```annotation_directory``` is also specified, predicted annotations supplement the existing annotations.

**Equivalent shell command**:
```
python panaroo_clean_inputs-runner.py -a <annotation_directory> -g <assembly_directory> -p <prodigal_directory> --index-file <index_file> -o {output} --threads <threads>
```

4. ```run_panaroo``` rule:

Runs panaroo on prokka-formatted functional annotation files. 

**Inputs**:
* ```input_directory```: Directory of prokka formatted function annotation files.
* ```output```: Directory for the panaroo outputs. 
* ```threads```: Number of threads for running panaroo. 

**Equivalent shell command**:
```
panaroo -i <input_directory>/*.gff -o <output> --clean-mode sensitive -t <threads>
```

 5. ```extract_assembly_stats``` rule:

Generates a JSON file of isolate metadata. 

**Inputs**:
* ```assembly_stats_directory```: Directory of uncompressed assembly statistics downloaded from NCBI GenBank.
* ```assembly_directory```: Directory of uncompressed assemblies downloaded from NCBI GenBank.
* ```index_file```: Filepath of JSON storing index numbers (Of the form ```{"isolateIndexNo": 0, "geneIndexNo": 0, "predictedIndexNo": 0}```).
* ```output```: Output directory for isolate metadata files. 
* ```email```: An email address to specify for Entrez programmatic access. 
* ```threads```: Number of threads for converting asssembly statistics to JSON.
* ```previous_run```: Directory storing previous snakemake outputs.  

**Outputs**:
*```isolateAssemblyAttributes.json```: JSON file of isolate metadata. 
*```biosampleIsolatePairs.json```: JSON file of isolate name key and BioSample accession ID values. 
*```indexIsolatePairs.json```: JSON file of isolate name key and isolate index number values. 

**Equivalent shell command**:
```
python extract_assembly_stats-runner.py -a <assembly_stats_directory> -g <assembly_directory> -i <index_file> -o <output> -e <email> --previous-run previous_run --threads <threads>```

