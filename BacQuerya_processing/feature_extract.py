#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses GFF and Sequence parsing to convert GFF files to json strings and extract gene names and sequences for COBS indexing.
"""
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from joblib import Parallel, delayed
import json
import glob
import gzip
import os
import sys
from tqdm import tqdm

def get_options():

    import argparse

    description = 'Extract features from gff and sequence files'
    parser = argparse.ArgumentParser(description=description,
                                        prog='feature_extract')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-s",
                        "--sequences",
                        dest="sequences",
                        required=True,
                        help='directory of genomic sequences',
                        type=str)
    io_opts.add_argument("-g",
                        "--gffs",
                        dest="gffs",
                        required=True,
                        help="directory of annotations in gff format",
                        type=str)
    io_opts.add_argument("-j",
                        "--isolate-json",
                        dest="isolate_json",
                        required=True,
                        help="json of all isolate attributes output by isolate_attributes",
                        type=str)
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_file",
                        required=True,
                        help="output file for gene json",
                        type=str)
    io_opts.add_argument("-i",
                        "--index-no",
                        dest="index_no",
                        required=False,
                        help="integer value to start gene index from",
                        default=0,
                        type=int)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    args = parser.parse_args()
    return (args)

def GFF_to_JSON(gff_file, seq_dir):
    """Use BCBio and SeqIO to convert isolate GFF files to JSON strings and identify the corresponding genomic sequence"""
    feature_list = []
    in_seq_file = os.path.basename(gff_file.replace(".gff.gz", ".fna.gz"))
    in_seq_file = os.path.join(seq_dir, in_seq_file)
    in_seq_handle = gzip.open(in_seq_file,'rt')
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()

    in_handle = gzip.open(gff_file,'rt')

    for rec in GFF.parse(in_handle, base_dict=seq_dict):
        sequence_record = rec.seq
        features = rec.features
        for f in features:
            start = int(f.location.start)
            end = int(f.location.end)
            strand = int(f.location.strand)
            if strand == 1:
                feature_sequence = str(sequence_record[start : end])
            elif strand == -1:
                feature_sequence = str(sequence_record[start : end].reverse_complement())
            json_features = {"type":f.type,
                            "location": {"strand":strand,
                                        "start":start,
                                        "end":end},
                            "sequence":feature_sequence,
                            "sequenceLength":len(feature_sequence)}
            qualifiers = f.qualifiers
            json_features.update(qualifiers)
            feature_list.append(json_features)
    in_handle.close()
    isolate_name = os.path.basename(gff_file).replace(".gff.gz", "")
    feature_dict = {"isolateName" : isolate_name.replace("_", " "),
                    "features" : feature_list}
    return feature_dict

def append_gene_indices(isolate_file, all_features):
    """Add indices of all genes within the isolate to the isolate attribute json"""
    with open(isolate_file, "r") as f:
        isolate_json = f.read()
    isolate_dict = json.loads(isolate_json)
    isolate_list_index = []
    for isol_name in isolate_dict["information"]:
        isolate_list_index.append(isol_name["isolateName"])
    for isol_features in all_features:
        isolate_gene_indices = []
        isolate_name = isol_features["isolateName"]
        for feature in isol_features["features"]:
            if feature["gbkey"][0] == "Gene" or feature["gbkey"][0] == "gene":
                isolate_gene_indices.append(feature["gene_index"])
        isolate_dict_position = isolate_list_index.index(isolate_name)
        isolate_dict["information"][isolate_dict_position]["geneIndices"] = isolate_gene_indices
    with open(isolate_file, "w") as o:
        o.write(json.dumps(isolate_dict))

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    output_dir = os.path.dirname(args.output_file)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    gffs = glob.glob(args.gffs + '/*.gff.gz')

    gff_list = [
        gffs[i:i + args.n_cpu] for i in range(0, len(gffs), args.n_cpu)
    ]
    # parrallelise feature extraction
    all_features = []
    all_genes = []
    for gff in tqdm(gff_list):
        features = Parallel(n_jobs=args.n_cpu)(delayed(GFF_to_JSON)(g,
                                                                    args.sequences) for g in gff)
        all_features += features
    index_no = args.index_no
    for single_isolate in all_features:
        feature_dict_list = single_isolate["features"]
        for json_features in feature_dict_list:
            if json_features["gbkey"][0] == "Gene" or json_features["gbkey"][0] == "gene":
                gene_dict = {"gene":json_features["Name"][0],
                            "sequence":json_features["sequence"],
                            "index":index_no}
                json_features["gene_index"] = index_no
                index_no += 1
                all_genes.append(gene_dict)
    # add gene indices to isolate jsons
    append_gene_indices(args.isolate_json,
                        all_features)
    with open(os.path.join(output_dir, "allIsolates.json"), "w") as a:
        a.write(json.dumps({"information":all_features}))
    with open(args.output_file, "w") as n:
        n.write(json.dumps({"information":all_genes}))
    sys.exit(0)

if __name__ == '__main__':
    main()
