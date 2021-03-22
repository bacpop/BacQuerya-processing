#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses GFF and Sequence parsing to convert GFF files to json strings and extract gene names and sequences for COBS indexing.
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
import pandas as pd
import re
import sys
from tqdm import tqdm

from BacQuerya_processing.index_gene_features import elasticsearch_isolates
from BacQuerya_processing.panaroo_clean_inputs import reverse_complement

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
                     isolateIndexJSON,
                     threads):
    """Extract all newly annotated/identified genes from panaroo graph and update geneJSON"""
    G = nx.read_gml(os.path.join(graph_dir, "final_graph.gml"))
    gene_data = pd.read_csv(os.path.join(graph_dir, "gene_data.csv"))
    # convert gene_data df to json to speed up sequence extraction
    sys.stderr.write('\nConverting gene data dataframe to JSON\n')
    gene_data_json = {}
    for row in tqdm(range(len(gene_data["clustering_id"]))):
        sequence = gene_data["dna_sequence"][row]
        cluster_dict = {gene_data["clustering_id"][row] : (gene_data["annotation_id"][row], sequence)}
        gene_data_json.update(cluster_dict)
    num_isolates = len(G.graph["isolateNames"])
    updated_genes = []
    panaroo_pairs = dict()
    name_set = set()
    sys.stderr.write('\nExtracting node information from Panaroo graph\n')
    for node in tqdm(G._node):
        y = G._node[node]
        frequency = round((len(y["members"])/num_isolates)*100, 1)
        member_labels = []
        gene_data_sequences = []
        annotation_ids = []
        for mem in range(len(y["members"])):
            member_labels.append(G.graph["isolateNames"][mem])
            gene_data_row = gene_data_json[y["geneIDs"].split(";")[mem]]
            member_node_sequence = gene_data_row[1]
            member_annotation_id = gene_data_row[0]
            gene_data_sequences.append(member_node_sequence)
            annotation_ids.append(member_annotation_id)
        isolate_indices = [isolateIndexJSON[label] for label in member_labels]
        panaroo_pairs.update({y["name"] : index_no})
        if not y["description"] == "":
            name_set.add(y["name"])
            gene_names = y["name"]
            updated_genes.append({"panarooNames" : gene_names,
                                  "panarooDescriptions" : y["description"].split(";"),
                                  "panarooFrequency": frequency,
                                  "gene_index": index_no,
                                  "foundIn_labels": member_labels,
                                  "foundIn_sequences": gene_data_sequences,
                                  "foundIn_indices": isolate_indices,
                                  "member_annotation_ids": annotation_ids})
        else:
            name_set.add(y["name"])
            gene_names = y["name"]
            updated_genes.append({"panarooNames" : gene_names,
                                  "panarooDescriptions" : ["hypothetical protein"],
                                  "panarooFrequency": frequency,
                                  "gene_index": index_no,
                                  "foundIn_labels": member_labels,
                                  "foundIn_sequences": gene_data_sequences,
                                  "foundIn_indices": isolate_indices,
                                  "member_annotation_ids": annotation_ids})
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
    json_feature_dicts = {}
    json_feature_list = []
    non_CDS_list = []
    with open(gff_file, "r") as g:
        gff_content = g.read().split("##")[2:-1]
    gff_region_title = ""
    for region in gff_content:
        if "sequence-region" == region.split(" ")[0]:
            gff_region_title = region.split(" ")[1]
        region_annotations = region.splitlines()
        for annotation_line in region_annotations:
            if gff_region_title == annotation_line.split("\t")[0]:
                region_sequence = seq_dict[gff_region_title]
                annotation_content = annotation_line.split("\t")
                qualifiers = annotation_content[8]
                start = int(annotation_content[3])
                end = int(annotation_content[4])
                strand = annotation_content[6]
                if strand == "+":
                    feature_sequence = str(region_sequence[start : end])
                elif strand == "-":
                    feature_sequence = str(region_sequence[start : end].reverse_complement())
                isolate_index = isolateIndexJSON[label]
                json_features = {"type":annotation_content[2],
                                "strand":strand,
                                "start":start,
                                "end":end,
                                "sequence":feature_sequence,
                                "sequenceLength":len(feature_sequence),
                                "isolateName": label,
                                "isolateIndex": isolate_index}
                qualifiers = qualifiers.split(";")
                for qual in qualifiers:
                    split = qual.split("=")
                    attribute = split[0]
                    val = split[1]
                    json_features.update({attribute : val})
                    json_feature_list.append(json_features)
    json_feature_dicts.update({label : json_feature_list})
    return json_feature_dicts

def append_gene_indices(isolate_file, all_features):
    """Add indices of all genes within the isolate to the isolate attribute json"""
    with open(isolate_file, "r") as f:
        isolate_json = f.read()
    isolate_dict = json.loads(isolate_json)
    for isol_name in tqdm(range(len(isolate_dict["information"]))):
        isolate_gene_indices = []
        isolate_gene_names = []
        isolate_non_CDS = []
        isolateMetadataName = isolate_dict["information"][isol_name]["isolateName"]
        for annotation_line in all_features:
            if annotation_line["isolateName"].replace("_", " ") == isolateMetadataName and "panarooNames" in annotation_line.keys():
                isolate_gene_indices.append(annotation_line["gene_index"])
                isolate_gene_names.append(annotation_line["panarooNames"])
            if annotation_line["isolateName"].replace("_", " ") == isolateMetadataName and "featureIndex" in annotation_line.keys():
                isolate_non_CDS.append(annotation_line["featureIndex"])
        if not len(isolate_gene_names) == 0:
            isolate_dict["information"][isol_name]["geneIndices"] = isolate_gene_indices
            isolate_dict["information"][isol_name]["panarooNames"] = isolate_gene_names
            isolate_dict["information"][isol_name]["nonCDSIndices"] = isolate_non_CDS
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
                                                         isolateIndexJSON,
                                                         args.n_cpu)
    else:
        updated_annotations = False
    gffs = glob.glob(args.gffs + '/*.gff')
    sys.stderr.write('\nConverting annotation files to JSON\n')
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
    sys.stderr.write('\nUpdating annotation JSON with Panaroo-sourced information\n')
    all_features = []
    for isolate_label, annotations in tqdm(json_feature_dicts.items()):
        non_CDS_index_no = index_no
        for annotation_line in annotations:
            if "ID" in annotation_line.keys() and "CDS" == annotation_line["type"]:
                for node in updated_annotations:
                    if annotation_line["ID"] in node["member_annotation_ids"]:
                        annotation_line.update(node)
                        all_features.append(annotation_line)
            else:
                annotation_line.update({"featureIndex": str(annotation_line["isolateIndex"]) + "_" + str(non_CDS_index_no)})
                all_features.append(annotation_line)
                non_CDS_index_no += 1
    # add gene indices to isolate jsons
    sys.stderr.write('\nAdding gene indices to isolate assembly JSONs\n')
    append_gene_indices(args.isolate_json,
                        all_features)
    if not args.elastic:
        sys.stderr.write('\nWriting gene JSON files\n')
        with open(os.path.join(args.output_dir, "annotatedNodes.json"), "w") as n:
            n.write(json.dumps({"information":all_features}))
    else:
        sys.stderr.write('\nBuilding Elastic Search index\n')
        elasticsearch_isolates(all_features,
                               args.index_name)
    sys.exit(0)

if __name__ == '__main__':
    main()
