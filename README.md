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
* ```threads```: Number of threads for retrieval. Entrez allows up to 3 queries/second without an API key and 10 queries/second with an API key. To specify an API key, see **Adding secrets**.

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

Runs panaroo [https://github.com/gtonkinhill/panaroo](https://github.com/gtonkinhill/panaroo) on prokka-formatted functional annotation files. 

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
