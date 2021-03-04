#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses Biopython ENTREZ to retrieve and download information of interest available through NCBI.
"""
from Bio import Entrez
from joblib import Parallel, delayed
import json
import os
import sys
import subprocess
from tqdm import tqdm
import urllib.request

def get_options():

    import argparse

    description = 'Download reads of interest by accession ID'
    parser = argparse.ArgumentParser(description=description,
                                        prog='extract_entrez_reads')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-s",
                        "--search_term",
                        dest="accessions",
                        required=True,
                        help='file of accessions for metadata download. One line per entry.',
                        type=str)
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,
                        help="output directory for file of interest",
                        type=str)
    io_opts.add_argument("-e",
                        "--email",
                        dest="email",
                        required=True,
                        help="specify email for entrez access",
                        type=str)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for retrieval",
                        default=1,
                        type=int)
    io_opts.add_argument("--number",
                        dest="number",
                        required=False,
                        help="maximum number of files to retrieve",
                        default=10000,
                        type=int)
    args = parser.parse_args()
    return (args)


def download_metadata(cleaned_accession,
                      email,
                      number,
                      output_dir):
    """Download the attribute for accessions of interest using Biopython Entrez"""
    failed_accessions = []
    Entrez.email = email
    handle = Entrez.read(Entrez.esearch(db="sra", term=cleaned_accession, retmax = number))
    assembly_ids = handle['IdList']
    try:
        esummary_handle = Entrez.esummary(db="biosample", id=assembly_ids[0], report="full")
        esummary_record = Entrez.read(esummary_handle, validate = False)
        accession_metadata = {cleaned_accession : dict(esummary_record["DocumentSummarySet"]["DocumentSummary"][0])}
        with open(os.path.join(output_dir, cleaned_accession + ".json"), "w") as f:
            f.write(json.dumps(accession_metadata))
    except ValueError:
        # the requested attribute is not present
        sys.stderr.write("Metadata is not available for: " + cleaned_accession)
        pass
    except:
        # issue with the request. Re-requesting often solves
        failed_accessions.append(cleaned_accession)
        sys.stderr.write("Request failed, the following accession will be re-requested: " + cleaned_accession)
    return failed_accessions

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    with open(args.accessions, "r") as s:
        sample_accessions = s.read()
    sample_list = sample_accessions.split("\n")
    cleaned_accessions = []
    for access in range(len(sample_list)):
        if not sample_list[access] == "":
            cleaned_accessions.append(sample_list[access])
    job_list = [
        cleaned_accessions[i:i + args.n_cpu] for i in range(0, len(cleaned_accessions), args.n_cpu)
    ]
    for job in tqdm(job_list):
        failed_accessions = Parallel(n_jobs=args.n_cpu)(delayed(download_metadata)(access,
                                                                                   args.email,
                                                                                   args.number,
                                                                                   args.output_dir) for access in job)
    failed_accessions = [failed for row in failed_accessions for failed in row]
    # ensure available attributes are downloaded for all accessions of interest
    for failed_access in failed_accessions:
        failed = download_metadata(failed_access,
                                   args.output_dir)
        if not len(failed) == 0:
            failed_accessions.append(failed)
    sys.exit(0)

if __name__ == '__main__':
    main()
