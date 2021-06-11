import glob
from joblib import Parallel, delayed
import json
import os
import ssl
import subprocess
import sys
from urllib.request import urlopen
from shutil import copyfileobj
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
    run:
        shell('python extract_entrez_information-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}')

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
    run:
        shell('python extract_entrez_information-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}')

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
    run:
        shell('python extract_entrez_information-runner.py -s {input} -e {params.email} --threads {params.threads} -o {output} -a {params.attribute}')

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
    params:
        threads=config['n_cpu']
    output:
        directory("prodigal_predicted_annotations")
    run:
        if not os.path.exists("prodigal_predicted_annotations2"):
            def multithread_prodigal(assembly, output_dir):
                output_file = os.path.join(output_dir, os.path.splitext(os.path.basename(assembly))[0] + ".gff")
                shell_command = "mkdir -p " + output_dir + " && prodigal -f gff -q -i " + assembly + " -o " + output_file
                subprocess.run(shell_command, shell=True, check=True)

            assemblies = glob.glob(os.path.join(input.genome_dir[0], "*.fna"))
            job_list = [
                assemblies[i:i + params.threads] for i in range(0, len(assemblies), params.threads)
            ]
            for job in tqdm(job_list):
                Parallel(n_jobs=params.threads)(delayed(multithread_prodigal)(assem,
                                                                            output[0]) for assem in job)
        else:
            shell("mkdir {output} && cp prodigal_predicted_annotations2/* {output}")

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
        threads=config["n_cpu"],
        run_type=config["run_type"]
    output:
        directory("panaroo_output")
    run:
        if os.path.exists("panaroo_output2"):
            shell("mkdir {output} && cp panaroo_output2/* {output}")
        else:
            if params.run_type == "reference":
                num_annotations = len(glob.glob(os.path.join(input[0], "*.gff")))
                if num_annotations > 1:
                    shell("panaroo -i {input}/*.gff -o {output} --clean-mode sensitive -t {params.threads}")
                else:
                    shell("mkdir {output}")
            else:
                shell("mkdir {output}")

# merge current panaroo output with previous panaroo outputs
rule merge_panaroo:
    input:
        current_output=rules.run_panaroo.output,
        annotation_dir=rules.reformat_annotations.output
    params:
        threads=config['n_cpu'],
        run_type=config["run_type"]
    output:
        touch("merge_panaroo.done")
    run:
        annotation_files = glob.glob(os.path.join(input.annotation_dir[0], "*.gff"))
        if params.run_type == "reference":
            if len(annotation_files) > 1:
                if os.path.exists("previous_run"):
                    shell("mkdir merged_panaroo_output && panaroo-merge -d {input.current_output} previous_run/panaroo_output -o merged_panaroo_output -t {params.threads} && cp -rf merged_panaroo_output/* panaroo_output && rm -rf merged_panaroo_output")
            if len(annotation_files) == 1:
                shell_command = " mkdir merged_panaroo_output && panaroo-integrate -d {input.current_output}/ -i " + annotation_files[0] + " -t {params.threads} -o merged_panaroo_output && && cp -rf merged_panaroo_output/* panaroo_output && rm -rf merged_panaroo_output"
                shell(shell_command)

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
        email=config['extract_entrez_information']['email'],
    run:
        shell('python extract_assembly_stats-runner.py -a {input.entrez_stats} -g {input.genome_files} -i {params.index} -o {output} -e {params.email} --previous-run previous_run --threads {params.threads}')


# retrieve ena read metadata
rule retrieve_ena_read_metadata:
    input:
        access_file=config['extract_read_metadata']['accession_file'],
        entrez_isolates=rules.extract_assembly_stats.output
    output:
        output_dir=directory('retrieved_ena_read_metadata'),
        run_accessions="retrieved_ena_read_metadata/fastq_links.txt",
        isolateJSON="retrieved_ena_read_metadata/isolateReadAttributes.json"
    params:
        index=config['extract_assembly_stats']['index_file'],
        email=config['extract_entrez_information']['email'],
        threads=config['n_cpu'],
        GPS=False
    run:
        if os.path.exists("retrieved_ena_read_metadata2"):
            shell("mkdir {output} && cp retrieved_ena_read_metadata2/* {output}")
        else:
            shell('python extract_read_metadata-runner.py -s {input.access_file} -r ena -i {params.index} -e {params.email} --previous-run previous_run --threads {params.threads} -o {output.output_dir}')

# retrieve raw reads from ENA
rule retrieve_ena_reads:
    input:
        rules.retrieve_ena_read_metadata.output.run_accessions
    params:
        threads=config['n_cpu']
    output:
        directory("retrieved_ena_reads")
    run:
        if os.path.exists("retrieved_ena_reads2"):
            shell("mkdir {output} && cp retrieved_ena_reads2/* {output}")
        else:
            def download_read(accession, output_dir):
                if "contigs" in accession:
                    ## currently only downloading assemblies and not read sets for efficiency
                    ssl._create_default_https_context = ssl._create_unverified_context
                    with urlopen(accession) as in_stream, open(os.path.join(output_dir, os.path.basename(accession)), 'wb') as out_file:
                        copyfileobj(in_stream, out_file)
                else:
                    pass
                return "success"

            os.mkdir(output[0])
            with open(input[0], "r") as f:
                run_accessions = f.read().splitlines()
            job_list = [
                run_accessions[i:i + params.threads] for i in range(0, len(run_accessions), params.threads)
            ]
            for job in tqdm(job_list):
                    results = Parallel(n_jobs=int(params.threads))(delayed(download_read)(access,
                                                                                          output[0]) for access in job)
        #'ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR164/ERR164407/ERR164407.fastq.gz {output}' ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR214/001/ERR2144781/ERR2144781_1.fastq.gz

# run mash screenon assemblies and reads
rule run_mash_screen:
    input:
        read_dir=rules.retrieve_ena_reads.output,
        assembly_dir=rules.unzip_genomes.output
    params:
        threads=config['n_cpu']
    output:
        directory("mash_results")
    run:
        def run_mash(sequenceFile, output_dir):
            """Multithreaded mash on genomic sequences"""
            shell_command = "./mash screen -w -p 3 refseq.genomes.k21s1000.msh " + sequenceFile + " > " + os.path.join(output_dir, os.path.splitext(os.path.basename(sequenceFile))[0]) + ".tab"
            subprocess.run(shell_command, shell=True, check=True)

        reads = glob.glob(os.path.join(input.read_dir[0], "*"))
        assemblies = glob.glob(os.path.join(input.assembly_dir[0], "*"))
        genomes = assemblies + reads
        os.mkdir(output[0])
        job_list = [
                genomes[i:i + params.threads] for i in range(0, len(genomes), params.threads)
            ]
        for job in tqdm(job_list):
            Parallel(n_jobs=params.threads)(delayed(run_mash)(sequence,
                                                              output[0]) for sequence in job)

# supplement isolate metadata with the mash screen output
rule supplement_isolate_metadata:
    input:
        mash_output=rules.run_mash_screen.output,
        ena_read_metadata=rules.retrieve_ena_read_metadata.output.output_dir,
        assembly_metadata=rules.extract_assembly_stats.output
    params:
        threads=config['n_cpu']
    output:
        touch("supplement_metadata.done")
    run:
        def appendMashAssemblies(isolate, mash_dir):
            """Add mash output to assembly metadata"""
            accession = isolate["isolateNameUnderscore"]
            with open(os.path.join(mash_dir, accession + ".tab")) as mashFile:
                mashOut = mashFile.read().splitlines()
            mashOut.sort(key=lambda x: x.split("\t")[0], reverse=True)
            mashIdentity = []
            mashHashes = []
            mashSpecies = []
            for row in mashOut:
                row = row.split("\t")
                # filter out all rows with fewer than 3 matching hashes and those of phages
                if int(row[1].split("/")[0]) > 2 and not "phage" in row[5]:
                        mashIdentity.append(row[0])
                        mashHashes.append(row[1])
                        mashSpecies.append(row[5])
            isolate["mashIdentity"] = mashIdentity
            isolate["mashHashes"] = mashHashes
            isolate["mashSpecies"] = mashSpecies
            return isolate

        def appendMashReads(read, mash_dir):
            """Add mash output to assembly metadata"""
            try:
                with open(os.path.join(mash_dir, read["BioSample"] + ".tab")) as mashFile:
                    mashOut = mashFile.read().splitlines()
            except:
                try:
                    with open(os.path.join(mash_dir, read["read_accession"] + ".tab")) as mashFile:
                        mashOut = mashFile.read().splitlines()
                except:
                    with open(os.path.join(mash_dir, read["run_accession"] + ".tab")) as mashFile:
                        mashOut = mashFile.read().splitlines()
            mashOut.sort(key=lambda x: x.split("\t")[0], reverse=True)
            mashIdentity = []
            mashHashes = []
            mashSpecies = []
            for row in mashOut:
                row = row.split("\t")
                # filter out all rows with fewer than 3 matching hashes and those of phages
                if int(row[1].split("/")[0]) > 2 and not "phage" in row[5]:
                    mashIdentity.append(row[0])
                    mashHashes.append(row[1])
                    mashSpecies.append(row[5])
            read["mashIdentity"] = mashIdentity
            read["mashHashes"] = mashHashes
            read["mashSpecies"] = mashSpecies
            return read

        with open(os.path.join(input.assembly_metadata[0],"isolateAssemblyAttributes.json"), "r") as assemblyFile:
            assembly_metadata = json.loads(assemblyFile.read())["information"]
        # iterate through isolates and import match output
        sys.stderr.write("\nAdding mash output to assemblies\n")
        job_list = [
                assembly_metadata[i:i + params.threads] for i in range(0, len(assembly_metadata), params.threads)
            ]
        updated_assembly_metadata = []
        for job in tqdm(job_list):
            updated_assembly_metadata += Parallel(n_jobs=params.threads)(delayed(appendMashAssemblies)(isolate,
                                                                                                       input.mash_output[0]) for isolate in job)
        with open(os.path.join(input.assembly_metadata[0], "isolateAssemblyAttributes.json"), "w") as assemblyOut:
            assemblyOut.write(json.dumps({"information": updated_assembly_metadata}))
        # iterate through isolates with reads
        with open(os.path.join(input.ena_read_metadata, "isolateReadAttributes.json"), "r") as readFile:
            ena_metadata = json.loads(readFile.read())["information"]
        sys.stderr.write("\nAdding mash output to reads\n")
        job_list = [
                ena_metadata[i:i + params.threads] for i in range(0, len(ena_metadata), params.threads)
            ]
        updated_ena_metadata = []
        for job in tqdm(job_list):
            updated_ena_metadata += Parallel(n_jobs=params.threads)(delayed(appendMashReads)(read,
                                                                                             input.mash_output[0]) for read in job)
        with open(os.path.join(input.ena_read_metadata, "isolateReadAttributes.json"), "w") as readOut:
            readOut.write(json.dumps({"information": updated_ena_metadata}))

# build gene JSONS from GFF and sequence files
rule extract_genes:
    input:
        annotations=rules.reformat_annotations.output,
        genomes=rules.unzip_genomes.output,
        assemblyStatDir=rules.extract_assembly_stats.output,
        graphDir=rules.run_panaroo.output,
        merged_panaroo=rules.merge_panaroo.output,
        mash_metadata=rules.supplement_isolate_metadata.output
    output:
        directory("extracted_genes")
    params:
        index=config['extract_assembly_stats']['index_file'],
        threads=config['n_cpu'],
        index_name=config['index_sequences']['elasticSearchIndex'],
        run_type=config["run_type"]
    run:
        if os.path.exists("extracted_genes2"):
            shell("mkdir {output} && cp extracted_genes2/* {output}")
        else:
            shell('python extract_genes-runner.py -s {input.genomes} -a {input.annotations} -m {input.assemblyStatDir} -g {input.graphDir} -i {params.index} -o {output} --threads {params.threads} --index-name {params.index_name} --prev-dir previous_run --run-type {params.run_type} --update')

# generate mafft alignments for panaroo output
rule mafft_align:
    input:
        merged_output=rules.merge_panaroo.output,
        extracted_genes=rules.extract_genes.output,
        graphDir=rules.run_panaroo.output
    params:
        threads=config['n_cpu']
    output:
        directory("aligned_gene_sequences")
    run:
        shell("python generate_alignments-runner.py --graph-dir {input.graphDir} --extracted-genes {input.extracted_genes} --output-dir {output} --threads {params.threads}")

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
        aligned_genes=rules.mafft_align.output,
        graphDir=rules.run_panaroo.output
    params:
        threads=config['n_cpu']
    output:
        touch("merge_runs.done")
    run:
        shell('python merge_runs-runner.py --ncbi-metadata {input.ncbiAssemblyStatDir} --graph-dir {input.graphDir} --geneMetadataDir {input.extractedGeneMetadata} --alignment-dir {input.aligned_genes} --accessionFile {input.currentRunAccessions} --previous-run previous_run --threads {params.threads}')

# build COBS index of gene sequences from the output of extract_genes
rule index_gene_sequences:
    input:
        merged_runs=rules.merge_runs.output,
        indexed_isolates=rules.index_isolate_attributes.output,
    output:
        directory("indexed_genes")
    params:
        k_mer=config['index_sequences']['kmer_length'],
        threads=config['n_cpu'],
        index_type=config['index_sequences']['gene_type'],
        elasticIndex=config['index_sequences']['elasticSearchIndex'],
        index_genes=config["index_genes"]
    run:
        if params.index_genes == True:
            shell('python index_gene_features-runner.py -t {params.index_type} -i previous_run/extracted_genes -g previous_run/panaroo_output -o {output} --kmer-length {params.k_mer} --threads {params.threads} --elastic-index --index {params.elasticIndex}')
        if params.index_genes == False:
            shell('python index_gene_features-runner.py -t {params.index_type} -i previous_run/extracted_genes -g previous_run/panaroo_output -o {output} --kmer-length {params.k_mer} --threads {params.threads}')

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
        indexed_isolates=rules.index_isolate_attributes.output,
        ena_metadata=rules.retrieve_ena_read_metadata.output.output_dir,
        mash_dir=rules.run_mash_screen.output,
        supplemented_metadata=rules.supplement_isolate_metadata.output
    shell:
        'rm -rf {input.supplemented_metadata} {input.mash_dir} {input.ena_metadata} {input.retrieved_genomes} {input.indexed_isolates} {input.prodigal_output} {input.aligned_genes} {input.retrieved_assembly_stats} {input.unzippedAnnotations} {input.retrieved_annotations} {input.unzipped_genomes} {input.mergeRuns} {input.panarooOutput} {input.extractedGeneMetadata} {input.ncbiAssemblyStatDir} {input.reformattedAnnotations} {input.merged_panaroo}'