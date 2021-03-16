#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses GFF and Sequence parsing to convert GFF files to json strings and extract gene names and sequences for COBS indexing.
Cannot currently deal with refound genes identified by Panaroo or unnamed genes.
"""
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from elasticsearch import Elasticsearch
from joblib import Parallel, delayed
import json
import glob
import networkx as nx
import os
import sys
from tqdm import tqdm

from BacQuerya_processing.index_gene_features import elasticsearch_isolates

def get_options():

    import argparse

    description = 'Extract features from gff and sequence files'
    parser = argparse.ArgumentParser(description=description,
                                        prog='extract_genes')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-s",
                        "--sequences",
                        dest="sequences",
                        required=True,
                        help='directory of genomic sequences',
                        type=str)
    io_opts.add_argument("-a",
                        "--gffs",
                        dest="gffs",
                        required=True,
                        help="directory of annotations in gff format",
                        type=str)
    io_opts.add_argument("-g",
                        "--graph-dir",
                        dest="graph_dir",
                        required=False,
                        default=False,
                        help="directory of Panaroo graph",
                        type=str)
    io_opts.add_argument("-k",
                        "--isolate-pairs",
                        dest="isolate_pairs",
                        required=True,
                        help="json file of index isolate name pairs",
                        type=str)
    io_opts.add_argument("-j",
                        "--isolate-json",
                        dest="isolate_json",
                        required=True,
                        help="json of all isolate attributes output by extract_assembly_stats",
                        type=str)
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,
                        help="output directory",
                        type=str)
    io_opts.add_argument("-i",
                        "--index-no",
                        dest="index_no",
                        required=False,
                        help="integer value to start gene index from",
                        default=0,
                        type=int)
    io_opts.add_argument("--elastic-index",
                        dest="elastic",
                        help="don't write gene json and index directly in script",
                        action='store_true',
                        default=False)
    io_opts.add_argument("--index-name",
                        dest="index_name",
                        required=False,
                        help="index to create/append to",
                        type=str)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    args = parser.parse_args()
    return (args)

def generate_library(graph_dir,
                     index_no,
                     output_dir,
                     isolateIndexJSON):
    """Extract all newly annotated/identified genes from panaroo graph and update geneJSON"""
    G = nx.read_gml(os.path.join(graph_dir, "final_graph.gml"))
    num_isolates = len(G.graph["isolateNames"])
    updated_genes = []
    panaroo_pairs = dict()
    name_set = set()
    for node in tqdm(G._node):
        y = G._node[node]
        frequency = round((len(y["members"])/num_isolates)*100, 1)
        member_labels = []
        for mem in y["members"]:
            member_labels.append(G.graph["isolateNames"][mem])
        panaroo_pairs.update({y["name"] : index_no})
        # if there are multiple sequ
        sequences = y["dna"]
        if not y["description"] == "":
            name_set.add(y["name"])
            gene_names = y["name"]
            updated_genes.append({"geneName" : gene_names,
                                  "description" : y["description"].split(";"),
                                  "geneFrequency": frequency,
                                  "gene_index": index_no,
                                  "members": member_labels})
        else:
            name_set.add(y["name"])
            gene_names = y["name"]
            updated_genes.append({"geneName" : gene_names,
                                  "description" : ["hypothetical protein"],
                                  "geneFrequency": frequency,
                                  "gene_index": index_no,
                                  "members": member_labels})
        index_no += 1
    # write name, index pairs in graph for COBS indexing in index_gene_features
    with open(os.path.join(output_dir, "panarooPairs.json"), "w") as o:
        o.write(json.dumps(panaroo_pairs))
    return updated_genes, index_no

def build_gff_jsons(gff_file,
                    seq_dir,
                    isolateIndexJSON):
    """Use BCBio and SeqIO to convert isolate GFF files to JSON strings and identify the corresponding genomic sequence"""
    label = os.path.basename(gff_file.replace(".gff", ""))
    in_seq_file = os.path.join(seq_dir, label + ".fna")
    in_seq_handle = open(in_seq_file,'r')
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()
    in_handle = open(gff_file,'r')
    json_feature_dicts = {}
    json_feature_list = []
    for rec in GFF.parse(in_handle, base_dict=seq_dict):
        sequence_record = rec.seq
        features = rec.features
        for f in features:
            qualifiers = f.qualifiers
            if "gene" in qualifiers.keys():
                start = int(f.location.start)
                end = int(f.location.end)
                strand = int(f.location.strand)
                if strand == 1:
                    feature_sequence = str(sequence_record[start : end])
                elif strand == -1:
                    feature_sequence = str(sequence_record[start : end].reverse_complement())
                isolate_index = isolateIndexJSON[label]
                json_features = {"type":f.type,
                                "strand":strand,
                                "start":start,
                                "end":end,
                                "sequence":feature_sequence,
                                "sequenceLength":len(feature_sequence),
                                "isolateName": label,
                                "isolateIndex": isolate_index}
                json_features.update(qualifiers)
                json_feature_list.append(json_features)
    json_feature_dicts.update({label : json_feature_list})
    in_handle.close()
    return json_feature_dicts

def build_panaroo_geneJSON(node_list, isolateIndexJSON, json_feature_dicts):
    "Build single layer JSON of genes and metadata in Panaroo graph for elasticsearch indexing"
    node_dicts = []
    for node in tqdm(node_list):
        isolate_labels = node["members"]
        isolate_indices = [isolateIndexJSON[label] for label in isolate_labels]
        nodeName = node["geneName"]
        gene_index = node["gene_index"]
        isolate_sequences = []
        isolates_annotated = []
        for isol in isolate_labels:
            if isol in json_feature_dicts.keys():
                annotations = json_feature_dicts[isol]
                for annot in annotations:
                    if "gene" in annot.keys() and annot["gene"][0] in nodeName.split("~~~"):
                        isolate_sequences.append(annot["sequence"])
                        isolates_annotated.append(annot["isolateName"])
        panaroo_dict = {"panarooNames": nodeName,
                        "panarooDescriptions": node["description"],
                        "panarooFrequency": node["geneFrequency"],
                        "gene_index": gene_index,
                        "foundIn_indices": isolate_indices,
                        "foundIn_labels": isolate_labels,
                        "annotatedIn_sequences": isolate_sequences,
                        "annotatedIn_labels": isolates_annotated}
        node_dicts.append(panaroo_dict)
    return node_dicts

def build_single_geneJSON():
    "Build single layer JSON of features not found in Panaroo graph for elasticsearch indexing"
    return

def append_gene_indices(isolate_file, all_features):
    """Add indices of all genes within the isolate to the isolate attribute json"""
    with open(isolate_file, "r") as f:
        isolate_json = f.read()
    isolate_dict = json.loads(isolate_json)
    for isol_name in tqdm(range(len(isolate_dict["information"]))):
        isolate_gene_indices = []
        isolate_gene_names = []
        isolateMetadataName = isolate_dict["information"][isol_name]["isolateName"]
        for annotation_line in all_features:
            if annotation_line["isolateName"].replace("_", " ") == isolateMetadataName and "panarooNames" in annotation_line.keys():
                isolate_gene_indices.append(annotation_line["gene_index"])
                isolate_gene_names.append(annotation_line["panarooNames"])
        if not len(isolate_gene_names) == 0:
            isolate_dict["information"][isol_name]["geneIndices"] = isolate_gene_indices
            isolate_dict["information"][isol_name]["panarooNames"] = isolate_gene_names
    with open(isolate_file, "w") as o:
        o.write(json.dumps(isolate_dict))

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # load isolate-index k, v pairs
    with open(args.isolate_pairs, "r") as isolateIndex:
        isolateString = isolateIndex.read()
    isolateIndexJSON = json.loads(isolateString)
    # standardise annotations from panaroo output
    index_no = args.index_no
    if args.graph_dir:
        sys.stderr.write('\nLoading Panaroo graph\n')
        updated_annotations, index_no = generate_library(args.graph_dir,
                                                         index_no,
                                                         args.output_dir,
                                                         isolateIndexJSON)
    else:
        updated_annotations = False
    gffs = glob.glob(args.gffs + '/*.gff')
    sys.stderr.write('\nConverting annotation files to JSONs\n')
    gff_list = [
        gffs[i:i + args.n_cpu] for i in range(0, len(gffs), args.n_cpu)
    ]
    # parrallelise conversion of GFF to json
    all_genes = []
    json_feature_dicts = {}
    for gff in tqdm(gff_list):
        features = Parallel(n_jobs=args.n_cpu)(delayed(build_gff_jsons)(g,
                                                                        args.sequences,
                                                                        isolateIndexJSON) for g in gff)
        for feature_item in features:
            json_feature_dicts.update(feature_item)
    sys.stderr.write('\nBuilding gene JSONs from Panaroo and annotation dicts\n')
    node_dicts = build_panaroo_geneJSON(updated_annotations,
                                        isolateIndexJSON,
                                        json_feature_dicts)
    # currently needed to index all genes when panaroo output is not available- will change to only features not in panaroo graph
    # iterate through isolate annotations to add the gene index and panaroo names
    sys.stderr.write('\nUpdating annotation JSON with Panaroo-sourced information\n')
    all_features = []
    for ndict in tqdm(node_dicts):
        panarooNames = ndict["panarooNames"]
        gene_index = ndict["gene_index"]
        for foundIn in ndict["foundIn_labels"]:
            foundIn_isolateAnnot = json_feature_dicts[foundIn]
            for annotation_line in foundIn_isolateAnnot:
                if "Name" in annotation_line.keys() and annotation_line["Name"][0] in panarooNames:
                    annotation_line.update({"panarooNames": panarooNames, "gene_index": gene_index})
                    all_features.append(annotation_line)
    # add gene indices to isolate jsons
    sys.stderr.write('\nAdding gene indices to isolate assembly JSONs\n')
    append_gene_indices(args.isolate_json,
                        all_features)
    if not args.elastic:
        sys.stderr.write('\nWriting gene JSON files\n')
        with open(os.path.join(args.output_dir, "annotatedNodes.json"), "w") as n:
            n.write(json.dumps({"information":node_dicts}))
    else:
        sys.stderr.write('\nBuilding Elastic Search index\n')
        elasticsearch_isolates(node_dicts,
                               args.index_name)
    sys.exit(0)

if __name__ == '__main__':
    main()
