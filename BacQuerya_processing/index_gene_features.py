"""
Build a searcheable COBS index from the output of extract_genes
- cannot index genes without indices because they are not in panaroo graph
"""
import cobs_index as cobs
from elasticsearch import Elasticsearch
import glob
from joblib import Parallel, delayed
import json
import networkx as nx
import os
import re
import shutil
import sys
from tqdm import tqdm
import tempfile

from BacQuerya_processing.secrets import ELASTIC_API_URL, ELASTIC_API_ID, ELASTIC_API_KEY

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
                        help="don't write gene json and index directly in script",
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
                           index_name):
    # rate of indexing decreases substantially after about 1500 items
    partioned_items = [
        allIsolatesJson[i:i + 1500] for i in range(0, len(allIsolatesJson), 1500)
        ]
    sys.stderr.write('\nIndexing CDS features\n')
    for item in tqdm(partioned_items):
        client = Elasticsearch([ELASTIC_API_URL],
                                api_key=(ELASTIC_API_ID, ELASTIC_API_KEY))
        # iterate through features
        for line in tqdm(item):
            if "gene_index" in line.keys():
                client = Elasticsearch([ELASTIC_API_URL],
                            api_key=(ELASTIC_API_ID, ELASTIC_API_KEY))
                response = client.index(index = index_name,
                                        id = line["gene_index"],
                                        body = line,
                                        request_timeout=30)
        #if "featureIndex" in line.keys():
           # response = client.index(index = index_name,
                                  #  id = line["featureIndex"],
                                  #  body = line,
                                  #  request_timeout=30)

def write_gene_files(gene_dict, temp_dir):
    """Write gene sequences to individual files with index as filename"""
    gene_index = str(gene_dict["gene_index"])
    gene_sequence = gene_dict["sequence"]
    with open(os.path.join(temp_dir, gene_index + ".txt"), "w") as g:
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
            with open(os.path.join(args.input_dir, "allIsolates.json"), "r") as inFeatures:
                geneString = inFeatures.read()
            isolateGeneDicts = json.loads(geneString)["information"]
            sys.stderr.write('\nBuilding elasticsearch index\n')
            elasticsearch_isolates(isolateGeneDicts,
                                args.index_name)
        # load panaroo graph and write sequence files from COG representatives
        sys.stderr.write('\nLoading panaroo graph\n')
        G = nx.read_gml(os.path.join(args.graph_dir, "final_graph.gml"))
        with open(os.path.join(os.path.dirname(args.graph_dir), "extracted_genes", "panarooPairs.json"), "r") as jsonFile:
            pairString = jsonFile.read()
        pairs = json.loads(pairString)
        panarooPairsUpdated = []
        representative_sequences = []
        sys.stderr.write('\nWriting gene-specific files for COBS indexing\n')
        for node in tqdm(G._node):
            y = G._node[node]
            gene_name = y["name"]
            dna = y["dna"].split(";")
            gene_index = pairs[y["name"]]
            for seq in range(len(dna)):
                representative_sequences.append({"gene_index": str(gene_index) + "_v" + str(seq), "sequence": dna[seq]})
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
