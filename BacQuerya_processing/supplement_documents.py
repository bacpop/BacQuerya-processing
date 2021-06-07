"""Supplement an entry in the elastic indices with the information from a JSON document. This document must have the gene or isolate index number as the key."""
from elasticsearch import Elasticsearch
from joblib import Parallel, delayed
import json
import sys
from tqdm import tqdm

from BacQuerya_processing.secrets import ELASTIC_API_URL, ELASTIC_ISOLATE_API_ID, ELASTIC_ISOLATE_API_KEY, ELASTIC_GENE_API_ID, ELASTIC_GENE_API_KEY

def get_options():

    import argparse
    description = 'Append information from a JSON document to information in an elastic index'
    parser = argparse.ArgumentParser(description=description,
                                     prog='supplement-entry')
    io_opts = parser.add_argument_group('Inputs')
    io_opts.add_argument("-i",
                        "--input",
                        dest="input_file",
                        required=False,
                        help='input JSON document from which to source the new information. The keys must nbe index numbers corresponding to documents in the index.',
                        type=str)
    io_opts.add_argument("--index",
                        dest="index_name",
                        required=False,
                        help="elasticsearch index to add to",
                        type=str,
                        choices=['gene_index', 'isolate_index'])
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    args = parser.parse_args()
    return (args)

def get_existing_entries(index_no, index_name, data):
    """Query elastic search index by document Id"""
    if index_name == "isolate":
        apiID = ELASTIC_ISOLATE_API_ID
        apiKEY = ELASTIC_ISOLATE_API_KEY
    if index_name == "gene":
        apiID = ELASTIC_GENE_API_ID
        apiKEY = ELASTIC_GENE_API_KEY
    client = Elasticsearch([ELASTIC_API_URL],
                           api_key=(apiID, apiKEY))
    fetchData = data[index_no]
    response = client.update(index=index_name, doc_type="_doc", id=index_no, body=fetchData)

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()
    # import the supplementary information
    with open(args.input_file, "r") as inJSON:
        supplement = json.loads(inJSON.read())
    # need to extract information already in the indices to prevent them being overwritten.
    indices = []
    for key in supplement.keys():
        indices.append(key)
    job_list = [
            indices[i:i + args.n_cpu] for i in range(0, len(indices), args.n_cpu)
            ]
    searchResults = []
    for job in tqdm(job_list):
        searchResults += Parallel(n_jobs=args.n_cpu)(delayed(get_existing_entries)(index,
                                                                                   args.index_name,
                                                                                   supplement) for index in job)
    sys.exit(0)

if __name__ == '__main__':
    main()
