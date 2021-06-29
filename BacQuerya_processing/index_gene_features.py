"""
Build a searcheable COBS index from the output of extract_genes
"""
#import cobs_index as cobs
from elasticsearch import Elasticsearch
import glob
from joblib import Parallel, delayed
import json
import networkx as nx
import os
import pyodbc
import re
import shutil
import sys
from tqdm import tqdm
import tempfile

from BacQuerya_processing.secrets import ELASTIC_API_URL, ELASTIC_GENE_API_ID, ELASTIC_GENE_API_KEY, SQL_CONNECTION_STRING

def get_options():

    import argparse
    description = 'Generate a searcheable COBS index from assemblies or gene sequences'
    parser = argparse.ArgumentParser(description=description,
                                     prog='index-genes')
    io_opts = parser.add_argument_group('Inputs')

    io_opts.add_argument("-t",
                        "--type",
                        dest="type",
                        required=True,
                        help="index assembly files or gene sequences",
                        choices=['assembly', 'gene'],
                        type=str)
    io_opts.add_argument("-i",
                        "--input-dir",
                        dest="input_dir",
                        required=False,
                        help='directory output by extract_genes (required for type=gene)',
                        type=str)
    io_opts.add_argument("-g",
                        "--graph-dir",
                        dest="graph_dir",
                        required=False,
                        help='directory of graph output by panaroo (required for type=gene and --all-genes=False)',
                        type=str)
    io_opts.add_argument("-a",
                        "--assembly-dir",
                        dest="assembly_dir",
                        required=False,
                        help='directory of assembly sequences (required for type=assembly)',
                        type=str)
    io_opts.add_argument("--isolate-dir",
                        dest="isolate_dir",
                        required=False,
                        help='directory of isolate metadata output by extract_assembly_stats.py',
                        type=str)
    io_opts.add_argument("-o",
                        "--output-dir",
                        dest="output_dir",
                        required=True,
                        help="output directory for generated index",
                        type=str)
    io_opts.add_argument("--kmer-length",
                        dest="kmer_length",
                        help="specify kmer length for the constructed index (default = 15)",
                        default=15,
                        type=int)
    io_opts.add_argument("--false-positive",
                        dest="fpr",
                        help="false positive rate for index. Greater fpr means smaller index (default = 0.01).",
                        default=0.01,
                        type=float)
    io_opts.add_argument("--index",
                        dest="index_name",
                        required=False,
                        help="elasticsearch index to create/append to",
                        type=str)
    io_opts.add_argument("--elastic-index",
                        dest="elastic",
                        help="index gene json in elastic index",
                        action='store_true',
                        default=False)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    io_opts.add_argument("--dirty",
                        dest="dirty",
                        help="keep gene files used to build the index",
                        action='store_true',
                        default=False)

    args = parser.parse_args()
    return (args)

def elasticsearch_isolates(allIsolatesJson,
                           index_name,
                           isolateMetadataDict):
    """Function to index gene metadata. Gene annotations are index using elastic and a list of isolates containing each gene is stored in a redis database"""
    # rate of indexing with elastic decreases substantially after about 1500 items
    partioned_items = [
        list(allIsolatesJson.keys())[i:i + 1500] for i in range(0, len(allIsolatesJson.keys()), 1500)
        ]
    sys.stderr.write('\nIndexing CDS features\n')
    with pyodbc.connect(SQL_CONNECTION_STRING) as conn:
        with conn.cursor() as cursor:
            cursor.execute('''CREATE TABLE GENE_METADATA
                (GENE_ID INT PRIMARY KEY     NOT NULL,
                 METADATA           TEXT    NOT NULL);''')
            #cursor.execute("DROP TABLE GENE_METADATA;")
    for keys in tqdm(partioned_items):
        elastic_client = Elasticsearch([ELASTIC_API_URL],
                                api_key=(ELASTIC_GENE_API_ID, ELASTIC_GENE_API_KEY))
        # iterate through features
        with pyodbc.connect(SQL_CONNECTION_STRING) as conn:
            with conn.cursor() as cursor:
                for line in tqdm(keys):
                    isolate_labels = allIsolatesJson[line]["foundIn_labels"]
                    isolate_indices = allIsolatesJson[line]["foundIn_indices"]
                    isolate_biosamples = allIsolatesJson[line]["foundIn_biosamples"]
                    isolate_annotationIDs = allIsolatesJson[line]["member_annotation_ids"]
                    # delete isolates from the gene metadata to reduce elastic entry size
                    del allIsolatesJson[line]["foundIn_labels"]
                    del allIsolatesJson[line]["foundIn_indices"]
                    del allIsolatesJson[line]["foundIn_biosamples"]
                    del allIsolatesJson[line]["member_annotation_ids"]
                    # extract supplementary metadata for isoaltes identified to be containing the gene
                    metaList = []
                    for isolate_row in isolateMetadataDict:
                        if isolate_row["isolate_index"] in isolate_indices:
                            isolate_metadata = {"BioSample": isolate_row["BioSample"],
                                                "sequenceURL": isolate_row["sequenceURL"]}
                            if "contig_stats" in isolate_row:
                                isolate_metadata.update({"contig_stats": isolate_row["contig_stats"]})
                            if "scaffold_stats" in isolate_row:
                                isolate_metadata.update({"scaffold_stats": isolate_row["scaffold_stats"]})
                            if "In_Silico_St" in isolate_row:
                                isolate_metadata.update({"In_Silico_St": isolate_row["In_Silico_St"]})
                            if "GPSC" in isolate_row:
                                isolate_metadata.update({"contig_stats": isolate_row["GPSC"]})
                            if "Country" in isolate_row:
                                isolate_metadata.update({"contig_stats": isolate_row["Country"]})
                            if "Year" in isolate_row:
                                isolate_metadata.update({"contig_stats": isolate_row["Year"]})
                            metaList.append(isolate_metadata)
                    response = elastic_client.index(index = index_name,
                                                    id = int(line),
                                                    body = allIsolatesJson[line],
                                                    request_timeout=60)
                    # store a list of isolates containing the gene in the SQL db
                    MetadataJSON = json.dumps({"foundIn_labels": isolate_labels,
                                "foundIn_indices": isolate_indices,
                                "foundIn_biosamples": isolate_biosamples,
                                "member_annotation_ids": isolate_annotationIDs,
                                "isolateMetadata": metaList})
                    db_command = "INSERT INTO GENE_METADATA (GENE_ID,METADATA) \
                        VALUES (" + str(line) + ", '" + MetadataJSON + "')"
                    cursor.execute(db_command)
    sys.stderr.write('\nGene metadata was indexed successfully\n')

def write_gene_files(gene_dict, temp_dir):
    """Write gene sequences to individual files with index as filename"""
    document_name = str(gene_dict["consistentNames"])
    gene_sequence = gene_dict["sequence"]
    with open(os.path.join(temp_dir, document_name + ".txt"), "w") as g:
        g.write(gene_sequence)

def write_assembly_files(assembly_file, temp_dir):
    """Write assembly sequences to individual files"""
    with open(assembly_file, "r") as f:
        assembly_sequence = f.read()
    assembly_sequence = re.sub(r'[^ACTG]', '', assembly_sequence)
    assembly_basename = os.path.basename(assembly_file).replace(".fna", ".txt")
    with open(os.path.join(temp_dir, assembly_basename), "w") as o:
        o.write(assembly_sequence)

def create_index(temp_dir, output_dir, kmer_length, fpr):
    """"Create COBS index"""
    params = cobs.CompactIndexParameters()
    params.term_size = kmer_length
    params.clobber = True               # overwrite output and temporary files
    params.false_positive_rate = fpr    # higher false positive rate -> smaller index
    cobs.compact_construct(temp_dir,
                           os.path.join(output_dir,
                                str(kmer_length) + "_index.cobs_compact"),
                           index_params=params)

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
    if args.type == "gene":
        if args.elastic:
            with open(os.path.join(args.input_dir, "annotatedNodes.json"), "r") as inFeatures:
                geneString = inFeatures.read()
            # supplement the gene metadata with metadata from the isolate metadata
            with open(os.path.join(args.isolate_dir, "isolateAssemblyAttributes.json"), "r") as inIsolates:
                isolateMetadataDict = json.loads(inIsolates.read())["information"]
            isolateGeneDicts = json.loads(geneString)["information"]
            sys.stderr.write('\nBuilding elasticsearch index\n')
            elasticsearch_isolates(isolateGeneDicts,
                                   args.index_name,
                                   isolateMetadataDict)
        # load panaroo graph and write sequence files from COG representatives
        sys.stderr.write('\nLoading panaroo graph\n')
        G = nx.read_gml(os.path.join(args.graph_dir, "final_graph.gml"))
        with open(os.path.join(os.path.dirname(args.graph_dir), "extracted_genes", "panarooPairs.json"), "r") as jsonFile:
            pairString = jsonFile.read()
        pairs = json.loads(pairString)
        representative_sequences = []
        sys.stderr.write('\nWriting gene-specific files for COBS indexing\n')
        # need to apply same constraints as those in extract_genes.py
        for node in G._node:
            y = G._node[node]
            gene_name = y["name"]
            splitNames = gene_name.split("~~~")
            # need to make sure we're not indexing annotations that have been predicted by prodigal and are not supported by existing annotations
            if not all("PRED_" in name for name in splitNames):
                dna = y["dna"].split(";")
                for key in pairs.keys():
                    if pairs[key]["panarooNames"] == y["name"]:
                        document_name = pairs[key]["consistentNames"]
                for seq in range(len(dna)):
                    representative_sequences.append({"consistentNames": str(document_name) + "_v" + str(seq), "sequence": dna[seq]})
        job_list = [
            representative_sequences[i:i + args.n_cpu] for i in range(0, len(representative_sequences), args.n_cpu)
        ]
        # parrallelise writing of gene-specific files for indexing
        for job in tqdm(job_list):
            Parallel(n_jobs=args.n_cpu)(delayed(write_gene_files)(feature,
                                                                  temp_dir) for feature in job)
    if args.type == "assembly":
        assembly_files_compressed = glob.glob(os.path.join(args.assembly_dir, "*.fna"))
        job_list = [
            assembly_files_compressed[i:i + args.n_cpu] for i in range(0, len(assembly_files_compressed), args.n_cpu)
        ]
        # parrallelise writing of assembly files to temp_dir for indexing
        for job in tqdm(job_list):
            Parallel(n_jobs=args.n_cpu)(delayed(write_assembly_files)(a,
                                                                      temp_dir) for a in job)
    sys.stderr.write('\nBuilding COBS sequence index\n')
    create_index(temp_dir,
                 args.output_dir,
                 args.kmer_length,
                 args.fpr)
    if not args.dirty:
        shutil.rmtree(temp_dir)
    sys.exit(0)

if __name__ == '__main__':
    main()
