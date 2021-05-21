#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Append isolate assembly stats to locally hosted elasticsearch index
"""
from elasticsearch import Elasticsearch
import json
import glob
import os
import sys
from tqdm import tqdm

from BacQuerya_processing.secrets import ELASTIC_API_URL, ELASTIC_ISOLATE_API_ID, ELASTIC_ISOLATE_API_KEY

def get_options():

    import argparse

    description = 'Index isolate json'
    parser = argparse.ArgumentParser(description=description,
                                        prog='index_isolate_attributes')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-i",
                        "--index",
                        dest="index",
                        required=True,
                        help='index to create/append to',
                        type=str)
    io_opts.add_argument("-f",
                        "--featureJSON",
                        dest="featureJSON",
                        required=False,
                        help='json file of assembled isolate features (extract_assembly_stats.py)',
                        type=str)
    io_opts.add_argument("-e",
                        "--ena-metadata-dir",
                        dest="ena_dir",
                        required=False,
                        help='directory of ENA read metadata (extract_read_metadata.py -r ena)',
                        type=str)
    io_opts.add_argument("-g",
                        "--feature-file",
                        dest="feature_file",
                        required=True,
                        help='json of all genes (extract_genes.py)',
                        type=str)
    args = parser.parse_args()
    return (args)

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(args.feature_file):
        raise AttributeError("Extract genes output is missing!")
    client = Elasticsearch([ELASTIC_API_URL],
                        api_key=(ELASTIC_ISOLATE_API_ID, ELASTIC_ISOLATE_API_KEY))
    if client.ping():
        sys.stderr.write('\nConnected to ES client\n')
    else:
        sys.stderr.write('\nCould not connect to ES client!\n')
    doc_list = []
    # index NCBI assembly features
    if args.featureJSON:
        with open(args.featureJSON, "r") as f:
            isolateString = f.read()
        # convert the string to a dict object
        dict_doc = json.loads(isolateString)
        doc_list += dict_doc['information']
    # index ENA read features
    if args.ena_dir:
        ena_metadata_file = os.path.join(args.ena_dir, "isolateReadAttributes.json")
        with open(ena_metadata_file, "r") as f:
            enaString = f.read()
            enaJSON = json.loads(enaString)
            doc_list += enaJSON["information"]
    for line in tqdm(doc_list):
        try:
            response = client.index(index = args.index,
                                    id = line["isolate_index"],
                                    body = line)
        except:
            sys.stderr.write('\nIssue indexing isolate: ' + line['isolateName'] + '\n')
    sys.exit(0)
