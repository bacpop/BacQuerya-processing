
configfile: 'config.yml'

# extract assembly-stats
rule retrieve_assembly_stats:
    input:
        config['entrez_extract']['accession_file']
    output:
        directory('retreived_assemblies')
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
        directory('retreived_annotations')
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
        directory('retreived_annotations')
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
        'retrieved_features/isolateFeatures.json'
    params:
        index=config['isolate_attributes']['index_no'],
        threads=config['n_cpu']
    shell:
       'python isolate_attributes-runner.py -a {input} -i {params.index} -o {output} --threads {params.threads}'

# build gene JSONS from GFF and sequence files
rule feature_extract:
    input:
        rules.retrieve_genomes.output
    output:
        directory("retrieved_genes")
    params:
        index=config['feature_extract']['index_no'],
        threads=config['n_cpu']
    shell:
       'python isolate_attributes-runner.py -a {input} -i {params.index} -o {output} --threads {params.threads}'

# append isolate attributes to elasticsearch index
rule index_isolates:
    input:
        rules.isolate_attributes.output
    params:
        index=config['index_isolate_attributes']['index']
    shell:
       'python index_isolate_attributes-runner.py -f {input} -i {params.index}'
