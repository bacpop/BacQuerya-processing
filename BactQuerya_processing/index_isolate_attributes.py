#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Append isolate features to locally hosted elasticsearch index
"""
from elasticsearch import Elasticsearch
import json
import os
import sys

def get_options():

    import argparse

    description = 'Index isolate json'
    parser = argparse.ArgumentParser(description=description,
                                        prog='index_isolate_features')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-f",
                        "--featureJSON",
                        dest="featureJSON",
                        required=True,
                        help='json file of isolate features',
                        type=str)
    io_opts.add_argument("-i",
                        "--index",
                        dest="index",
                        required=True,
                        help='index to create/append to',
                        type=str)
    io_opts.add_argument("-g",
                        "--feature-file",
                        dest="feature_file",
                        required=True,
                        help='file of all genes output by feature_extract',
                        type=str)
    args = parser.parse_args()
    return (args)

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.isfile(args.feature_file):
        raise AttributeError("Feature file is missing!")
    client = Elasticsearch("localhost:9200")
    if client.ping():
        sys.stderr.write('\nConnected to ES client\n')
    else:
        sys.stderr.write('\nCould not connect to ES client!\n')

    with open(args.featureJSON, "r") as f:
        isolateString = f.read()
    # convert the string to a dict object
    dict_doc = json.loads(isolateString)
    doc_list = dict_doc['information']
    for line in doc_list:
        try:
            response = client.index(index = args.index,
                                    id = line["index"],
                                    body = line)
        except:
            sys.stderr.write('\nIssue indexing isolate: ' + line['isolateName'] + '\n')
    sys.exit(0)
