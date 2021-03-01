#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract attributes from isolate-specific assembly stats and constuct a JSON for all attributes in all isolates.
"""
from assembly_stats import read_genome, calculate_stats
from joblib import Parallel, delayed
import json
import glob
import os
import re
import sys
from tqdm import tqdm

def get_options():

    import argparse

    description = 'Extract assembly stats'
    parser = argparse.ArgumentParser(description=description,
                                        prog='extract_assembly_stats')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-a",
                        "--assemblies",
                        dest="assemblies",
                        required=True,
                        help='directory of assembly reports',
                        type=str)
    io_opts.add_argument("-g",
                        "--genomes",
                        dest="genomes",
                        required=True,
                        help='directory of genomic sequences',
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

def calculate_assembly_stats(genomeFile):
    contig_lens, scaffold_lens, gc_cont = read_genome(genomeFile)
    contig_stats = calculate_stats(contig_lens, gc_cont)
    scaffold_stats = calculate_stats(scaffold_lens, gc_cont)
    return contig_stats, scaffold_stats

def assembly_to_JSON(assigned_index, genome_dir):
    """Use assembly stats to extract information for elasticsearch indexing"""
    index_no = assigned_index['isolate_index']
    assembly_file = assigned_index['assembly file']
    isolate_name = os.path.basename(assembly_file).replace("_assembly_stats.txt", "")
    genome_file = os.path.join(genome_dir, isolate_name + ".fna")
    contig_stats, scaffold_stats = calculate_assembly_stats(genome_file)
    with open(assembly_file, "r") as f:
        assembly_features = f.read().split("\n")
    assembly_dict = {"isolateName" : isolate_name.replace("_", " "),
                    "isolate_index" : index_no,
                    "contig_stats": contig_stats,
                    "scaffold_stats": scaffold_stats}
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
        assigned_index = {"isolate_index": index_no, "assembly file": assembly}
        index_no += 1
        indexed_assemblies.append(assigned_index)

    job_list = [
        indexed_assemblies[i:i + args.n_cpu] for i in range(0, len(indexed_assemblies), args.n_cpu)
    ]
    # parrallelise assembly feature extraction
    all_features = []
    for job in tqdm(job_list):
        features = Parallel(n_jobs=args.n_cpu)(delayed(assembly_to_JSON)(assem, args.genomes) for assem in job)
        all_features += features
    with open(os.path.join(args.output_file), "w") as a:
        a.write(json.dumps({"information":all_features}))
    sys.exit(0)

if __name__ == '__main__':
    main()
