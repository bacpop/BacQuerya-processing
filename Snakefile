
configfile: 'config.yml'

# extract assembly-stats
rule retrieve_assembly_stats:
    input:
        config['entrez_extract']['accession_file']
    output:
        directory('retrieved_assemblies')
    params:
        email=config['entrez_extract']['email'],
        attribute=config['entrez_extract']['assembly'],
        threads=config['n_cpu']
    shell:
       'python entrez_extract-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}'

# extract isolate-GFFS
rule retrieve_annotations:
    input:
        config['entrez_extract']['accession_file']
    output:
        directory('retrieved_annotations')
    params:
        email=config['entrez_extract']['email'],
        attribute=config['entrez_extract']['gff'],
        threads=config['n_cpu']
    shell:
       'python entrez_extract-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}'

# extract isolate-genomes
rule retrieve_genomes:
    input:
        config['entrez_extract']['accession_file']
    output:
        directory('retrieved_genomes')
    params:
        email=config['entrez_extract']['email'],
        attribute=config['entrez_extract']['genome'],
        threads=config['n_cpu']
    shell:
       'python entrez_extract-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}'

# build isolate JSONS from assembly-stats
rule isolate_attributes:
    input:
        rules.retrieve_assembly_stats.output
    output:
        'isolate_attributes/isolateAttributes.json'
    params:
        index=config['isolate_attributes']['index_no'],
        threads=config['n_cpu']
    shell:
       'python isolate_attributes-runner.py -a {input} -i {params.index} -o {output} --threads {params.threads}'

# build gene JSONS from GFF and sequence files
rule feature_extract:
    input:
        annotations=rules.retrieve_annotations.output,
        genomes=rules.retrieve_genomes.output,
        isolateJson=rules.isolate_attributes.output,
    output:
        "isolate_genes/allGenes.json"
    params:
        index=config['feature_extract']['index_no'],
        threads=config['n_cpu']
    shell:
       'python feature_extract-runner.py -s {input.genomes} -g {input.annotations} -j {input.isolateJson} -i {params.index} -o {output} --threads {params.threads}'

# append isolate attributes to elasticsearch index
rule index_isolate_attributes:
    input:
        attribute_file=rules.isolate_attributes.output,
        feature_file=rules.feature_extract.output
    params:
        index=config['index_isolate_attributes']['index'],
    shell:
       'python index_isolate_attributes-runner.py -f {input.attribute_file} -i {params.index} -g {input.feature_file}'

# build COBS index of genetic sequences from the output of feature_extract
rule index_gene_features:
    input:
        rules.feature_extract.output
    output:
        directory("sequence_index")
    params:
        k_mer=config['index_gene_features']['kmer_length'],
        threads=config['n_cpu']
    shell:
       'python index_gene_features-runner.py -i {input} -o {output} --kmer-length {params.k_mer} --threads {params.threads}'
