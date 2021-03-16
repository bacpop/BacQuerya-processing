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

# retrieve raw reads using SRA toolkit
rule retrieve_sra_reads:
    input:
        config['extract_read_metadata']['accession_file']
    output:
        directory('retrieved_sra_reads')
    shell:
       'prefetch --output-directory {output} -v --option-file {input}'

# retrieve sra read metadata
rule retrieve_sra_read_metadata:
    input:
        config['extract_read_metadata']['accession_file']
    output:
        output_dir=directory('retrieved_sra_read_metadata'),
    params:
        email=config['extract_entrez_information']['email'],
        threads=config['n_cpu']
    shell:
       'python extract_read_metadata-runner.py -s {input} -r sra -e {params.email} --threads {params.threads} -o {output.output_dir}'

# retrieve ena read metadata
rule retrieve_ena_read_metadata:
    input:
        config['extract_read_metadata']['accession_file']
    output:
        output_dir=directory('retrieved_ena_read_metadata'),
        run_accessions="retrieved_ena_read_metadata/fastq_links.txt",
        isolateJSON="retrieved_ena_read_metadata/isolateReadAttributes.json"
    params:
        email=config['extract_entrez_information']['email'],
        threads=config['n_cpu']
    shell:
       'python extract_read_metadata-runner.py -s {input} -r ena -i 200 -e {params.email} --threads {params.threads} -o {output.output_dir}'

# retrieve raw reads from ENA
rule retrieve_ena_reads:
    input:
        rules.retrieve_ena_read_metadata.output.run_accessions
    output:
        directory("retrieved_ena_reads")
    run:
        with open(input[0], "r") as f:
            run_accessions = f.read().splitlines()
        for access in run_accessions:
            shell_command = "wget --directory-prefix " + output[0] + " " + access
            shell(shell_command)
    #'ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR164/ERR164407/ERR164407.fastq.gz {output}' ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR214/001/ERR2144781/ERR2144781_1.fastq.gz

# gunzip genome files
rule unzip_genomes:
    input:
        rules.retrieve_genomes.output
    output:
        directory("unzipped_genomes")
    shell:
        "mkdir {output} && cp {input}/*.gz {output} && gunzip {output}/*.gz"

# convert .sra to fastq
rule expand_sra_reads:
    input:
        rules.retrieve_sra_reads.output
    output:
        directory("sra_reads")
    shell:
        "fasterq-dump --split-files -O {output} {input}/*.sra"

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
            shell_command = "meryl k=" + best_kmer_result + " count output " + meryl_foldername + " " + assem
            shell(shell_command)

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
        shell_command = "$MERQURY/merqury.sh meryl_merged_files " + assembly_string + " output"
        shell(shell_command)

# clean up merqury outputs
rule clean_merqury_outputs:
    input:
        merqury_output=rules.run_merqury.output.output_dir
    shell:
        "mv output.* completeness.stats *.gz {input.merqury_output} && rm -rf *.gz*"

# run prodigal to predict genes in assemblies
rule run_prodigal:
    input:
        rules.unzip_genomes.output
    output:
        directory("prodigal_predicted_annotations")
    run:
        assemblies = glob.glob(os.path.join(input[0], "*.fna"))
        for assembly in assemblies:
            output_file = os.path.join(output[0], os.path.splitext(os.path.basename(assembly))[0])
            shell("mkdir -p {output} && prodigal -f gff -i " + assembly + " -o " + output_file + ".gff")

# reformat annotation files for panaroo input
rule reformat_annotations:
    input:
        genome_dir=rules.unzip_genomes.output,
        annotation_dir=rules.unzip_annotations.output
    params:
        threads=config['n_cpu']
    output:
        directory("panaroo_cleaned_annotations")
    shell:
        "python panaroo_clean_inputs-runner.py -a {input.annotation_dir} -g {input.genome_dir} -o {output} --threads {params.threads}"

# run panaroo on reformatted annotations
rule run_panaroo:
    input:
        rules.reformat_annotations.output
    params:
        threads=config['n_cpu']
    output:
        directory("panaroo_output")
    shell:
        "panaroo -i {input}/*.gff -o {output} --clean-mode sensitive -t {params.threads}"

# build isolate JSONS from assembly-stats
rule extract_assembly_stats:
    input:
        entrez_stats=rules.retrieve_assembly_stats.output,
        genome_files=rules.unzip_genomes.output
    output:
        isolateFile='extracted_assembly_stats/isolateAssemblyAttributes.json',
        indexJSON='extracted_assembly_stats/indexIsolatePairs.json'
    params:
        index=config['extract_assembly_stats']['index_no'],
        threads=config['n_cpu']
    shell:
       'python extract_assembly_stats-runner.py -a {input.entrez_stats} -g {input.genome_files} -i {params.index} -o {output.isolateFile} -k {output.indexJSON} --threads {params.threads}'

# build gene JSONS from GFF and sequence files
rule extract_genes:
    input:
        annotations=rules.unzip_annotations.output,
        genomes=rules.unzip_genomes.output,
        isolateJson=rules.extract_assembly_stats.output.isolateFile,
        #graphDir=rules.run_panaroo.output
        graphDir="panaroo_merged_output",
        isolateKeyPairs=rules.extract_assembly_stats.output.indexJSON
    output:
        directory("extracted_genes")
    params:
        index=config['extract_genes']['index_no'],
        threads=config['n_cpu'],
        index_name=config['index_sequences']['elasticSearchIndex']
    shell:
       'python extract_genes-runner.py -s {input.genomes} -a {input.annotations} -g {input.graphDir} -j {input.isolateJson} -k {input.isolateKeyPairs} -i {params.index} -o {output} --threads {params.threads} --elastic-index --index-name {params.index_name}'

# build gene JSONS from prodigal-predicted GFF and sequence files
rule extract_predicted_genes:
    input:
        annotations=rules.run_prodigal.output,
        #genomes=rules.assembled_reads.output,
        isolateJson=rules.retrieve_ena_read_metadata.output.isolateJSON,
    output:
        directory("extracted_genes")
    params:
        index=config['extract_genes']['index_no'],
        threads=config['n_cpu']
    shell:
       'python extract_genes-runner.py -s {input.genomes} -g {input.annotations} -j {input.isolateJson} -i {params.index} -o {output} --threads {params.threads}'

# append isolate attributes to elasticsearch index
rule index_isolate_attributes:
    input:
        attribute_file=rules.extract_assembly_stats.output,
        feature_file=rules.extract_genes.output,
        ena_metadata=rules.retrieve_ena_read_metadata.output.output_dir
    params:
        index=config['index_isolate_attributes']['index'],
    shell:
       'python index_isolate_attributes-runner.py -f {input.attribute_file} -e {input.ena_metadata} -i {params.index} -g {input.feature_file}/allGenes.json'

# build COBS index of gene sequences from the output of extract_genes
rule index_gene_sequences:
    input:
        input_dir=rules.extract_genes.output,
        graph_dir="panaroo_merged_output"
    output:
        directory("index_genes")
    params:
        k_mer=config['index_sequences']['kmer_length'],
        threads=config['n_cpu'],
        index_type=config['index_sequences']['gene_type'],
        elasticIndex=config['index_sequences']['elasticSearchIndex']
    shell:
       'python index_gene_features-runner.py -t {params.index_type} -i {input.input_dir} -g {input.graph_dir} -o {output} --kmer-length {params.k_mer} --threads {params.threads} --index {params.elasticIndex}'

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
       'python index_gene_features-runner.py -t {params.index_type} -a {input} -o {output} --kmer-length {params.k_mer} --threads {params.threads}'
