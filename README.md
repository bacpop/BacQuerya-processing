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

### 1. retrieve_assembly_stats, retrieve_annotations and retrieve_genomes rules:

Retrieves information of interest from NCBI GenBank by accession ID. 

**Inputs**:
* ```accession_file```: Filepath of a "\n" separated list of BioSample or assembly accession IDs for asssemblies available through NCBI GenBank or Refseq. 
* ```attribute```: Retrieve assembly sequences, functional annotations or assembly statistics (specified as ```genomes```, ```annotation``` or ```assembly-stats```).
* ```email```: An email address to specify for Entrez programmatic access. 
* ```output```: Output directory name for retrieved attributes. 
* ```threads```: Number of threads for retrieval. Entrez allows up to 3 queries/second without an API key and 10 queries/second with an API key. To specify an API key, see **Local instances of BacQuerya**.

**Equivalent shell command**:
```
python extract_entrez_information-runner.py -s <accession_file> -e <email> --threads <threads> -o <output> -a <attribute>
```

### 2. run_prodigal rule:

Predicts genes in all assemblies in an input directory and outputs predicted annotations in GFF3 format.
**Inputs**:
* ```genome_dir```: Directory containing assemblies. 
* ```output```: Output directory name for predicted annotations. 
* ```threads```: Number of threads for prediction. 

### 3. eformat_annotations rule:

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

### 4. run_panaroo rule:

Runs panaroo ([https://github.com/gtonkinhill/panaroo](https://github.com/gtonkinhill/panaroo)) on prokka-formatted functional annotation files. 

**Inputs**:
* ```input_directory```: Directory of prokka formatted function annotation files.
* ```output```: Directory for the panaroo outputs. 
* ```threads```: Number of threads for running panaroo. 

**Equivalent shell command**:
```
panaroo -i <input_directory>/*.gff -o <output> --clean-mode sensitive -t <threads>
```

 ### 5. extract_assembly_stats rule:

Generates a JSON file of metadata for isolates with assemblies. 

**Inputs**:
* ```assembly_stats_directory```: Directory of uncompressed assembly statistics downloaded from NCBI GenBank.
* ```assembly_directory```: Directory of uncompressed assemblies downloaded from NCBI GenBank.
* ```index_file```: Filepath of JSON storing index numbers (Of the form ```{"isolateIndexNo": 0, "geneIndexNo": 0, "predictedIndexNo": 0}```).
* ```output```: Output directory for isolate metadata files. 
* ```email```: An email address to specify for Entrez programmatic access. 
* ```threads```: Number of threads for converting asssembly statistics to JSON.
* ```previous_run```: Directory storing previous snakemake outputs.  

**Outputs**:
* ```isolateAssemblyAttributes.json```: JSON file of isolate assembly metadata. 
* ```biosampleIsolatePairs.json```: JSON file of isolate name key and BioSample accession ID values. 
* ```indexIsolatePairs.json```: JSON file of isolate name key and isolate index number values. 

**Equivalent shell command**:
```
python extract_assembly_stats-runner.py -a <assembly_stats_directory> -g <assembly_directory> -i <index_file> -o <output> -e <email> --previous-run previous_run --threads <threads>
```

### 6. retrieve_ena_read_metadata rule:

Generates a JSON file of metadata for isolates with reads. 

**Inputs**:
* ```accession_file```: Filepath of a "\n" separated list of ERR or ERS accession IDs for asssemblies available through the ENA.  
* ```index_file```: Filepath of JSON storing index numbers (Of the form ```{"isolateIndexNo": 0, "geneIndexNo": 0, "predictedIndexNo": 0}```).
* ```output```: Output directory for isolate metadata files. 
* ```email```: An email address to specify for Entrez programmatic access. 
* ```threads```: Number of threads for converting asssembly statistics to JSON.
* ```previous_run```: Directory storing previous snakemake outputs.  

**Outputs**:
* ```fastq_links.txt```: "\n" separated list of read sequence download URLs. 
* ```isolateReadAttributes.json```: JSON file of isolate metadata. 

**Equivalent shell command**:
```
python extract_read_metadata-runner.py -s <accession_file> -r ena -i <index_file> -o <output> -e <email> --previous-run previous_run --threads <threads>
```

### 7. retrieve_ena_reads rule:

Multithreaded download of read sets. 

**Inputs**:
* ```fastq_links```: Filepath of a "\n" separated list of read sequence download URLs. 

### 8. extract_genes rule:

Extracts gene metadata for assemblies from a panaroo output and adds genes per isolate to a JSON file of isolate metadata.

**Inputs**:
* ```assembly_directory```: Directory of uncompressed assemblies downloaded from NCBI GenBank.
* ```annotation_directory```: Directory of uncompressed functional annotations in GFF3 format.
* ```graph_directory```: Directory of a panaroo graph constructed from the annotations in the ```annotation_directory```.
* ```isolate_metadata```: Directory of isolate metadata output by ```extract_assembly_stats.py```
* ```index_metadata```: Directly index gene metadata in elastic index in this script (True or False).
* ```output```: Output directory for gene metadata files. 
* ```threads```: Number of threads to annotate "query isolates". 
* ```index_name```: Name of elastic gene metadata index.
* ```run_type```: Whether this is a "reference" or "query" run (see **reference vs query runs**). 
* ```update```: Update the input panaroo outputs with the new gene names and annotations. 

**Equivalent shell command**:
```
python extract_genes-runner.py -s <assembly_directory> -a <annotation_directory> -g <graph_directory> -m <isolate_metadata> -i <index_metadata> -o <output> --threads <threads> --index-name <index_name> --prev-dir previous_run --run-type <run_type> [--update]
```

### 9. mafft_align rule:

Extracts gene metadata for assemblies from a panaroo output and adds genes per isolate to a JSON file of isolate metadata.

**Inputs**:
* ```graph_directory```: Directory of a panaroo graph constructed from the annotations in the ```annotation_directory``` and updated by the ```extract_genes``` rule.
* ```gene_metadata```: Directory output by the ```extract_genes``` rule. Requires ```panarooPairs.json``` to ensure gene names are consistent across outputs.
* ```output```: Output directory for the aligned genes.
* ```threads```: Number of threads for alignment.

**Equivalent shell command**:
```
python generate_alignments-runner.py --graph-dir <graph_directory> --extracted-genes <gene_metadata> --output-dir <output> --threads <threads>
```

### 10. index_isolate_attributes rule:

Indexes isolate assembly and read metadata in a searchable elastic index. 

**Inputs**:
* ```isolate_asssembly_metadata```: Filepath of JSON of isolate assembly metadata.
* ```isolate_read_metadata```: Directory of isolate read metadata output by the ```retrieve_ena_read_metadata``` rule.
* ```index_name```: Name of elastic isolate metadata index.
* ```gene_metadata```: Directory output by the ```extract_genes``` rule. Required to ensure genes contained have been added to the isolate assembly metadata.

**Equivalent shell command**:
```
python index_isolate_attributes-runner.py -f <isolate_asssembly_metadata> -e <isolate_read_metadata> -i <index_name> -g <gene_metadata>
```

### 11. index_gene_sequences rule:

Indexes gene metadata in a searchable elastic index. 

**Inputs**:
* ```gene_metadata```: Directory containing extracted gene metadata. If using the snakemake pipeline, this will be ```previous_run/extracted_genes``` as the metadata for current and previous snakemake runs will have been merged. 
* ```graph_directory```: Directory of a panaroo graph constructed from the annotations in the ```annotation_directory``` and updated by the ```extract_genes``` rule. If using the snakemake pipeline, this will be ```previous_run/panaroo_output``` as the metadata for current and previous snakemake runs will have been merged. 
* ```output```: Output directory for the constructed COBS index.
* ```kmer_length```: K-mer length at which to construct the COBS index (default=31).
* ```threads```: Number of threads to write fasta files for COBS input.
* ```index_name```: Name of elastic isolate metadata index.
* ```elastic-index```: Index genes in elastic gene index. If excluded, only a COBS gene index is constructed.  

**Equivalent shell command**:
```
python index_gene_features-runner.py -t gene -i <gene_metadata> -g <graph_directory> -o <output> --kmer-length <kmer_length> --threads <threads> --index  <index_name> [--elastic-index] 
```

## reference vs query runs

Users can choose the ```run_type``` for the BacQuerya Snakemake pipeline as a ```reference``` or ```query``` run with:
```
snakemake --cores 1 --config run_type=<run_type>
```
This option impacts whether or not the genomic information in the current run is used to update the existing panaroo output of a previous ```reference``` run. If ```reference``` is specified, panaroo is run on the isolates in the current run and this output is merged with the previous output to become the new reference panaroo output. If ```query``` is selected, prokka-formatted functional annotations are integrated into the reference output one by one to identify genes, but the reference panaroo output is **not** updated with the genomic information from these isolates. 

# Local instances of BacQuerya

BacQuerya-processing primarily has been designed to populate the indices for the BacQuerya website but we appreciate our processing and indexing pipeline may have use cases beyond our original intentions. This may include locally hosting indices, allowing users to index and search through sensitive genomic data not currently in the public domain. To do this, we anticipate that a number of scripts will require some customisation depending on the use case and we would be happy to answer any questions you may have on setting up one of these private instances. However, local instances of BacQuerya sourcing publicly available data can be set up relatively easily.

## Locally hosting elastic indices

Elasticsearch is easy to install, free for local installations and interacts with the processing pipeline through the same API as the elastic cloud indices we use with BacQuerya. Elasticsearch is available for download from [https://www.elastic.co/downloads/elasticsearch](https://www.elastic.co/downloads/elasticsearch). Uncompress the download and add the following line to the end of the ```config/elasticsearch.yml``` file within the bundle. 
```
http.cors:
  enabled: true
  allow-origin: /https?:\/\/localhost(:[0-9]+)?/
```
To run Elasticsearch at [http://localhost:9200](http://localhost:9200), enter the bundle directory and start the Elastic instance by running:
* on windows
```
bin\elasticsearch.bat
```
* on macOS/Linux
```
bin/elasticsearch
```

# Specifying elastic parameters and secrets

Our Elastic parameters and API keys are not available for public use therefore, to make BacQuerya-processing communicate with your local elasticsearch instance you must create a ```secrets.py``` file in the ```BacQuerya_processing``` directory and define the following paramters for export:

* ```ELASTIC_API_URL```: The elasticsearch endpoint. This is [http://localhost:9200](http://localhost:9200) for local instances. 
* ```ELASTIC_ISOLATE_API_ID```: The API ID for an API key with write access for your elastic isolate index. 
* ```ELASTIC_ISOLATE_API_KEY```: The API key with write access for your elastic isolate index. 
* ```ELASTIC_GENE_API_ID```: The API ID for an API key with write access for your elastic gene index. 
* ```ELASTIC_GENE_API_KEY```: The API key with write access for your elastic gene index.

To create API keys for your indices, see [https://www.elastic.co/guide/en/elasticsearch/reference/current/security-api-create-api-key.html](https://www.elastic.co/guide/en/elasticsearch/reference/current/security-api-create-api-key.html). Users may also define an ```ENTREZ_API_KEY``` here to increase the number of possible Entrez queries per second.

# Contributors

BacQuerya, BacQuerya-api and BacQuerya-processing were developed by Daniel Anderson and John Lees.
