#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses GFF parsing to convert GFF files to json strings.
"""
from BCBio import GFF
from Bio import SeqIO
from joblib import Parallel, delayed
import json
import glob
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
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,
                        help="output directory for gene json",
                        type=str)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    args = parser.parse_args()
    return (args)

def GFF_to_JSON(gff_file, output_dir):
    """Use BCBio and SeqIO to convert isolate GFF files to JSON strings"""
    feature_list = []

    in_seq_file = gff_file.replace(".gff", ".fna")
    in_seq_handle = open(in_seq_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()

    in_handle = open(gff_file)

    for rec in GFF.parse(in_handle, base_dict=seq_dict):
        features = rec.features
        for f in features:
            json_features = {"type":f.type,
                            "location": {"strand":f.location.strand,
                                        "start":int(f.location.start),
                                        "end":int(f.location.end)},
                            "id":f.id,
                            "qualifiers":f.qualifiers}
            feature_list.append(json_features)
    in_handle.close()
    isolate_name = os.path.basename(gff_file).replace(".gff", "")
    feature_dict = {"isolateName" : isolate_name.replace("_", " "),
                    "features" : feature_list}
    with open(os.path.join(output_dir, isolate_name + ".json"), "w") as f:
        f.write(json.dumps(feature_dict))
    return feature_dict

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    gffs = glob.glob(args.gffs + '/*.gff')

    job_list = [
        gffs[i:i + args.n_cpu] for i in range(0, len(gffs), args.n_cpu)
    ]
    # parrallelise feature extraction
    all_features = []
    all_isolate_names = []
    for job in tqdm(job_list):
        features = Parallel(n_jobs=args.n_cpu)(delayed(GFF_to_JSON)(g,
                                                                    args.output_dir) for g in job)
        all_features += features
    for i in all_features:
        all_isolate_names.append(i["isolateName"])
    with open(os.path.join(args.output_dir, "allIsolates.json"), "w") as a:
        a.write(json.dumps({"information":all_features}))
    with open(os.path.join(args.output_dir, "isolateNames.json"), "w") as n:
        n.write(json.dumps({"information":all_isolate_names}))
    sys.exit(0)

if __name__ == '__main__':
    main()
