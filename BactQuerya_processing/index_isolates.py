#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract features from isolate-specific assembly stats and constuct a JSON for all features in all isolates.
"""
from joblib import Parallel, delayed
import json
import glob
import os
import re
import sys
from tqdm import tqdm

def get_options():

    import argparse

    description = 'Extract features from gff and sequence files'
    parser = argparse.ArgumentParser(description=description,
                                        prog='feature_extract')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-a",
                        "--assemblies",
                        dest="assemblies",
                        required=True,
                        help='directory of assembly reports',
                        type=str)
    io_opts.add_argument("-i",
                        "--index-no",
                        dest="index_no",
                        required=False,
                        help="integer value to start index from",
                        default=0,
                        type=int)
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_file",
                        required=True,
                        help="output file for json of isolates",
                        type=str)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    args = parser.parse_args()
    return (args)

def assembly_to_JSON(assigned_index):
    """Use assembly stats to extract information for elasticsearch indexing"""
    index_no = assigned_index['index']
    assembly_file = assigned_index['assembly file']
    with open(assembly_file, "r") as f:
        assembly_features = f.read().split("\n")
    isolate_name = os.path.basename(assembly_file).replace("_assembly_stats.txt", "")
    assembly_dict = {"isolateName" : isolate_name.replace("_", " "),
                    "index" : index_no}
    for line in assembly_features:
        try:
            attribute = re.search('# (.*?):', line).group(1).replace(" ", "_")
            feature_dict = {attribute : line.split(":")[-1].strip().replace("_", " ")}
            assembly_dict.update(feature_dict)
        except AttributeError:
            pass
    return assembly_dict

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(os.path.dirname(args.output_file)):
        os.mkdir((os.path.dirname(args.output_file)))
    assembly_reports = glob.glob(args.assemblies + '/*_assembly_stats.txt')

    indexed_assemblies = []
    index_no = args.index_no
    for assembly in assembly_reports:
        assigned_index = {"index": index_no, "assembly file": assembly}
        index_no += 1
        indexed_assemblies.append(assigned_index)

    job_list = [
        indexed_assemblies[i:i + args.n_cpu] for i in range(0, len(indexed_assemblies), args.n_cpu)
    ]
    # parrallelise assembly feature extraction
    all_features = []
    for job in tqdm(job_list):
        features = Parallel(n_jobs=args.n_cpu)(delayed(assembly_to_JSON)(assem) for assem in job)
        all_features += features
    with open(os.path.join(args.output_file), "w") as a:
        a.write(json.dumps({"information":all_features}))
    sys.exit(0)

if __name__ == '__main__':
    main()
