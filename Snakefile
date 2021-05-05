import glob
import os
import subprocess
from tqdm import tqdm

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
        index=config['extract_assembly_stats']['index_file'],
        email=config['extract_entrez_information']['email'],
        threads=config['n_cpu']
    shell:
       'python extract_read_metadata-runner.py -s {input} -r ena -i {params.index} -e {params.email} --threads {params.threads} -o {output.output_dir}'

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
        genome_dir=rules.unzip_genomes.output,
        annotation_dir=rules.unzip_annotations.output
    output:
        directory("prodigal_predicted_annotations")
    run:
        assemblies = glob.glob(os.path.join(input.genome_dir[0], "*.fna"))
        existing_annotations = glob.glob(os.path.join(input.annotation_dir[0], "*.gff"))
        existing_annotations = [os.path.basename(filename).split(".gff")[0] for filename in existing_annotations]
        for assembly in tqdm(assemblies):
            #if not os.path.basename(assembly).split(".fna")[0] in existing_annotations:
            output_file = os.path.join(output[0], os.path.splitext(os.path.basename(assembly))[0] + ".gff")
            shell("mkdir -p {output} && prodigal -f gff -q -i " + assembly + " -o " + output_file)

# reformat annotation files for panaroo input
rule reformat_annotations:
    input:
        genome_dir=rules.unzip_genomes.output,
        annotation_dir=rules.unzip_annotations.output,
        prodigal_dir=rules.run_prodigal.output,
        indexFile=config['extract_assembly_stats']['index_file']
    params:
        threads=config['n_cpu']
    output:
        directory("panaroo_cleaned_annotations")
    shell:
        "python panaroo_clean_inputs-runner.py -a {input.annotation_dir} -g {input.genome_dir} -p {input.prodigal_dir} --index-file {input.indexFile} -o {output} --threads {params.threads}"

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

# generate mafft alignments for panaroo output
rule mafft_align:
    input:
        rules.run_panaroo.output
    params:
        threads=config['n_cpu']
    output:
        directory("aligned_gene_sequences")
    shell:
        "python generate_alignments-runner.py --graph-dir {input} --output-dir {output} --threads {params.threads}"

# merge current panaroo output with previous panaroo outputs
rule merge_panaroo:
    input:
        current_output=rules.run_panaroo.output
    params:
        threads=config['n_cpu']
    output:
        touch("merge_panaroo.done")
    run:
        if os.path.exists("previous_run"):
            shell("panaroo-merge -d {input.current_output} previous_run/panaroo_output -o previous_run/panaroo_output -t {params.threads}")
        else:
            os.mkdir("previous_run")
            shell("cp -r {input.current_output} previous_run")

# build isolate JSONS from assembly-stats
rule extract_assembly_stats:
    input:
        entrez_stats=rules.retrieve_assembly_stats.output,
        genome_files=rules.unzip_genomes.output
    output:
        directory("extracted_assembly_stats")
    params:
        index=config['extract_assembly_stats']['index_file'],
        threads=config['n_cpu'],
        email=config['extract_entrez_information']['email']
    shell:
       'python extract_assembly_stats-runner.py -a {input.entrez_stats} -g {input.genome_files} -i {params.index} -o {output}/isolateAssemblyAttributes.json -k {output}/indexIsolatePairs.json -b {output}/biosampleIsolatePairs.json -e {params.email} --threads {params.threads}'

# build gene JSONS from GFF and sequence files
rule extract_genes:
    input:
        annotations=rules.unzip_annotations.output,
        genomes=rules.unzip_genomes.output,
        assemblyStatDir=rules.extract_assembly_stats.output,
        merged_panaroo=rules.merge_panaroo.output
    output:
        directory("extracted_genes")
    params:
        index=config['extract_assembly_stats']['index_file'],
        threads=config['n_cpu'],
        index_name=config['index_sequences']['elasticSearchIndex']
    shell:
       'python extract_genes-runner.py -s {input.genomes} -a {input.annotations} -g previous_run/panaroo_output -j {input.assemblyStatDir}/isolateAssemblyAttributes.json -k {input.assemblyStatDir}/indexIsolatePairs.json -i {params.index} -o {output} --threads {params.threads} --elastic-index --index-name {params.index_name} --biosampleJSON {input.assemblyStatDir}/biosampleIsolatePairs.json'

# append isolate attributes to elasticsearch index
rule index_isolate_attributes:
    input:
        assemblyStatDir=rules.extract_assembly_stats.output,
        feature_file=rules.extract_genes.output,
        ena_metadata=rules.retrieve_ena_read_metadata.output.output_dir
    params:
        index=config['index_isolate_attributes']['index'],
    output:
        touch("index_isolates.done")
    shell:
       'python index_isolate_attributes-runner.py -f {input.assemblyStatDir}/isolateAssemblyAttributes.json -e {input.ena_metadata} -i {params.index} -g {input.feature_file}'

# merge the current run with information from previous runs
rule merge_runs:
    input:
        ncbiAssemblyStatDir=rules.extract_assembly_stats.output,
        extractedGeneMetadata=rules.extract_genes.output,
        currentRunAccessions=config['extract_entrez_information']['accession_file'],
        aligned_genes=rules.mafft_align.output
    params:
        threads=config['n_cpu']
    output:
        touch("merge_runs.done")
    shell:
        'python merge_runs-runner.py --ncbi-metadata {input.ncbiAssemblyStatDir} --geneMetadataDir {input.extractedGeneMetadata} --alignment-dir {input.aligned_genes} --accessionFile {input.currentRunAccessions} --previous-run previous_run --threads {params.threads}'

# build COBS index of gene sequences from the output of extract_genes
rule index_gene_sequences:
    input:
        input_dir=rules.extract_genes.output,
        merged_runs=rules.merge_runs.output,
        merged_panaroo=rules.merge_panaroo.output,
        fake_input=rules.index_isolate_attributes.output
    output:
        directory("index_genes")
    params:
        k_mer=config['index_sequences']['kmer_length'],
        threads=config['n_cpu'],
        index_type=config['index_sequences']['gene_type'],
        elasticIndex=config['index_sequences']['elasticSearchIndex'],
    shell:
       'python index_gene_features-runner.py -t {params.index_type} -i {input.input_dir} -g previous_run/panaroo_output -o {output} --kmer-length {params.k_mer} --threads {params.threads}'

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

# run entire pipeline and delete current run output directories when done
rule run_pipeline:
    input:
        ncbiAssemblyStatDir=rules.extract_assembly_stats.output,
        extractedGeneMetadata=rules.extract_genes.output,
        panarooOutput=rules.run_panaroo.output,
        mergeRuns=rules.merge_runs.output,
        reformattedAnnotations=rules.reformat_annotations.output,
        unzippedAnnotations=rules.unzip_annotations.output,
        retrieved_annotations=rules.retrieve_annotations.output,
        unzipped_genomes=rules.unzip_genomes.output,
        retrieved_genomes=rules.retrieve_genomes.output,
        retrieved_assembly_stats=rules.retrieve_assembly_stats.output,
        merged_panaroo=rules.merge_panaroo.output,
        aligned_genes=rules.mafft_align.output,
        prodigal_output=rules.run_prodigal.output,
        indexed_gene_sequences=rules.index_gene_sequences.output,
        indexed_isolates=rules.index_isolate_attributes.output
    shell:
        'rm -rf {input.retrieved_genomes} {input.indexed_isolates} {input.prodigal_output} {input.aligned_genes} {input.retrieved_assembly_stats} {input.unzippedAnnotations} {input.retrieved_annotations} {input.unzipped_genomes} {input.mergeRuns} {input.panarooOutput} {input.extractedGeneMetadata} {input.ncbiAssemblyStatDir} {input.reformattedAnnotations} {input.merged_panaroo}'