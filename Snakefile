import glob
import os
import subprocess

configfile: 'config.yml'

# retrieve assembly-stats
rule retrieve_assembly_stats:
    input:
        config['extract_entrez_information']['accession_file']
    output:
        directory('retrieved_assemblies')
    params:
        email=config['extract_entrez_information']['email'],
        attribute=config['extract_entrez_information']['assembly'],
        threads=config['n_cpu']
    shell:
       'python extract_entrez_information-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}'

# retrieve isolate-GFFS
rule retrieve_annotations:
    input:
        config['extract_entrez_information']['accession_file']
    output:
        directory('retrieved_annotations')
    params:
        email=config['extract_entrez_information']['email'],
        attribute=config['extract_entrez_information']['gff'],
        threads=config['n_cpu']
    shell:
       'python extract_entrez_information-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}'

# gunzip annotation files
rule unzip_annotations:
    input:
        rules.retrieve_annotations.output
    output:
        directory("unzipped_annotations")
    shell:
        "mkdir {output} && cp {input}/*.gz {output} && gunzip {output}/*.gz"

# retrieve isolate-genomes
rule retrieve_genomes:
    input:
        config['extract_entrez_information']['accession_file']
    output:
        directory('retrieved_genomes')
    params:
        email=config['extract_entrez_information']['email'],
        attribute=config['extract_entrez_information']['genome'],
        threads=config['n_cpu']
    shell:
       'python extract_entrez_information-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}'

# gunzip genome files
rule unzip_genomes:
    input:
        rules.retrieve_genomes.output
    output:
        directory("unzipped_genomes")
    shell:
        "mkdir {output} && cp {input}/*.gz {output} && gunzip {output}/*.gz"

# build single meryl dbs
rule single_meryl_dbs:
    input:
        genomes=rules.unzip_genomes.output
    output:
        single_files=directory(config["single_meryl_dbs"]["single_files"]),
        assembly_txt=config["single_meryl_dbs"]["assembly_txt_file"]
    run:
        # extract assembly statistics for illumina whole genome sequencing reads
        assembly_files =  glob.glob(os.path.join(input.genomes[0], "*.gz"))
        single_meryl_dir = output.single_files
        if not os.path.exists(single_meryl_dir):
            os.mkdir(single_meryl_dir)
        # write out assembly string for run_merqury
        assembly_string = " ".join(assembly_files)
        with open(output.assembly_txt, "w") as a:
            a.write(assembly_string)

        for assem in assembly_files:
            with open(assem, "rt") as f:
                genome = f.read().splitlines()
            genome_len = 0
            for line in genome:
                if not ">" in line:
                    genome_len += len(line)

            # calculate best k-mer length for genome
            best_kmer_result = subprocess.check_output("sh $MERQURY/best_k.sh " + str(genome_len), shell=True)
            best_kmer_result = best_kmer_result.decode('utf-8')
            best_kmer_result = best_kmer_result.splitlines()[2]
            # build meryl db for each assembly
            meryl_foldername = os.path.join(single_meryl_dir, os.path.splitext(os.path.splitext(os.path.basename(assem))[0])[0] + ".meryl")
            subprocess.run("meryl k=" + best_kmer_result + " count output " + meryl_foldername + " " + assem, shell=True, check=True)

# merge single meryl dbs for merqury input
rule merge_single_meryl_dbs:
    input:
        rules.single_meryl_dbs.output.single_files
    output:
        directory(config["merge_single_meryl_dbs"]["merged_db_folder"])
    shell:
        'meryl union-sum output {output} {input}/*.meryl'

# evaluate assembly quality using merqury
rule run_merqury:
    input:
        merged_dbs=rules.merge_single_meryl_dbs.output,
        assembly_txt=rules.single_meryl_dbs.output.assembly_txt
    output:
        output_dir=directory(config["run_merqury"]["merqury_output"])
    run:
        if not os.path.exists(output.output_dir):
            os.mkdir(output.output_dir)
        with open(input.assembly_txt, "r") as f:
            assembly_string = f.read()
        subprocess.run("$MERQURY/merqury.sh meryl_merged_files " + assembly_string + " output", shell=True, check=True)

# clean up merqury outputs
rule clean_merqury_outputs:
    input:
        merqury_output=rules.run_merqury.output.output_dir
    shell:
        "mv output.* completeness.stats *.gz {input.merqury_output} && rm -rf *.gz*"

# build isolate JSONS from assembly-stats
rule extract_assembly_stats:
    input:
        entrez_stats=rules.retrieve_assembly_stats.output,
        genome_files=rules.unzip_genomes.output
    output:
        'extracted_assembly_stats/isolateAttributes.json'
    params:
        index=config['extract_assembly_stats']['index_no'],
        threads=config['n_cpu']
    shell:
       'python extract_assembly_stats-runner.py -a {input.entrez_stats} -g {input.genome_files} -i {params.index} -o {output} --threads {params.threads}'

# build gene JSONS from GFF and sequence files
rule extract_genes:
    input:
        annotations=rules.unzip_annotations.output,
        genomes=rules.unzip_genomes.output,
        isolateJson=rules.extract_assembly_stats.output,
    output:
        "extracted_genes/allGenes.json"
    params:
        index=config['extract_genes']['index_no'],
        threads=config['n_cpu']
    shell:
       'python extract_genes-runner.py -s {input.genomes} -g {input.annotations} -j {input.isolateJson} -i {params.index} -o {output} --threads {params.threads}'

# append isolate attributes to elasticsearch index
rule index_isolate_attributes:
    input:
        attribute_file=rules.extract_assembly_stats.output,
        feature_file=rules.extract_genes.output
    params:
        index=config['index_isolate_attributes']['index'],
    shell:
       'python index_isolate_attributes-runner.py -f {input.attribute_file} -i {params.index} -g {input.feature_file}'

# build COBS index of gene sequences from the output of extract_genes
rule index_gene_sequences:
    input:
        rules.extract_genes.output
    output:
        directory("index_genes")
    params:
        k_mer=config['index_sequences']['kmer_length'],
        threads=config['n_cpu'],
        index_type=config['index_sequences']['gene_type']
    shell:
       'python index_gene_features-runner.py -t {params.index_type} -f {input} -o {output} --kmer-length {params.k_mer} --threads {params.threads}'

# build COBS index of gene sequences from the output of extract_genes
rule index_assembly_sequences:
    input:
        rules.unzip_genomes.output
    output:
        directory("index_assemblies")
    params:
        k_mer=config['index_sequences']['kmer_length'],
        threads=config['n_cpu'],
        index_type=config['index_sequences']['assembly_type']
    shell:
       'python index_gene_features-runner.py -t {params.index_type} -d {input} -o {output} --kmer-length {params.k_mer} --threads {params.threads}'
