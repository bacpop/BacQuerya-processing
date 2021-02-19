
configfile: 'config.yml'

# extract assembly-stats
rule retrieve_assembly_stats:
    input:
        config['entrez_extract']['accession_file']
    output:
        directory('retreived_entries')
    params:
        email=config['entrez_extract']['email'],
        attribute=config['entrez_extract']['attribute'],
        threads=config['n_cpu']
    shell:
       'python entrez_extract-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}'

# build isolate JSONS from assembly-stats
rule isolate_features:
    input:
        rules.retrieve_assembly_stats.output
    output:
        'retrieved_features/isolateFeatures.json'
    params:
        index=config['index_isolates']['index_no'],
        threads=config['n_cpu']
    shell:
       'python isolate_features-runner.py -a {input} -i {params.index} -o {output} --threads {params.threads}'

# append isolate features to elasticsearch index
rule index_isolates:
    input:
        rules.isolate_features.output
    params:
        index=config['isolate_assembly_features']['index']
    shell:
       'python index_isolate_features-runner.py -f {input} -i {params.index}'
