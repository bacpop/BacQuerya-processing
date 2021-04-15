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

    description = 'Download information of interest from NCBI'
    parser = argparse.ArgumentParser(description=description,
                                        prog='entrez_extract')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-s",
                        "--search_term",
                        dest="accessions",
                        required=True,
                        help='file of sequences to download. One line per entry. Specify accessions of interest or "genus species" for all accessions.',
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
    io_opts.add_argument("-a",
                        "--attributes",
                        dest="attrs",
                        required=False,
                        help="specify attributes to download",
                        default="genome",
                        choices=['genome', 'annotation', 'assembly-report', 'assembly-stats'],
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

def translate_attr(attribute):
    """Determine attribute suffixes for download from NCBI"""
    attr_dict = {'genome': ['_genomic.fna.gz', '.fna.gz'],
                'annotation':['_genomic.gff.gz', '.gff.gz'],
                'assembly-report':['_assembly_report.txt', '_assembly_report.txt'],
                'assembly-stats':['_assembly_stats.txt', '_assembly_stats.txt']}
    entrez_attribute = attr_dict[attribute][0]
    attribute_suffix = attr_dict[attribute][1]
    return entrez_attribute, attribute_suffix

def download_entries(cleaned_accession,
                    entrez_attribute,
                    attribute_suffix,
                    email,
                    number):
    """Download the attribute for accessions of interest using Biopython Entrez"""
    successful_accessions = []
    Entrez.email = email
    handle = Entrez.read(Entrez.esearch(db="assembly", term=cleaned_accession, retmax = number))
    assembly_ids = handle['IdList']
    try:
        esummary_handle = Entrez.esummary(db="assembly", id=assembly_ids[0], report="full")
        esummary_record = Entrez.read(esummary_handle, validate = False)
        url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        # search in GenBank if no RefSeq entry found
        if url == '':
            url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        label = os.path.basename(url)
        if entrez_attribute == "_assembly_stats.txt":
            sequenceLinks = {"accession": cleaned_accession,
                             "sequenceURL": url.replace("ftp://", "https://") + "/" + label + "_genomic.fna.gz"}
            with open(label + "_additionalAssemblyStats.txt", "w") as a:
                a.write(json.dumps(sequenceLinks))
        attribute_link = os.path.join(url,label + entrez_attribute)
        urllib.request.urlretrieve(attribute_link, label + attribute_suffix)
        successful_accessions.append(cleaned_accession)
    except ValueError:
        # the requested attribute is not present
        sys.stderr.write("\nThe requested attribute is not present for: " + cleaned_accession + "\n")
        successful_accessions.append(cleaned_accession)
    except:
        # issue with the request. Re-requesting often solves
        sys.stderr.write("\nRequest failed, the following accession will be re-requested: " + cleaned_accession + "\n")
    return successful_accessions

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

    entrez_attribute, attribute_suffix = translate_attr(args.attrs)

    # set wd to output_dir for urllib
    os.chdir(args.output_dir)
    sys.stderr.write("\nDownloading information for " + str(len(cleaned_accessions)) + " isolates\n")
    job_list = [
        cleaned_accessions[i:i + args.n_cpu] for i in range(0, len(cleaned_accessions), args.n_cpu)
    ]
    successful_accessions = []
    for job in tqdm(job_list):
        successful_accessions += Parallel(n_jobs=args.n_cpu)(delayed(download_entries)(str(access),
                                                                                       entrez_attribute,
                                                                                       attribute_suffix,
                                                                                       args.email,
                                                                                       args.number) for access in job)
    # ensure available attributes are downloaded for all accessions of interest
    successful_accessions = set(successful_accessions)
    while len(successful_accessions) != len(cleaned_accessions):
        for access in cleaned_accessions:
            if not access in successful_accessions:
                sys.stderr.write("\nRerequesting isolate: " + str(failed_access) + "\n")
                success = download_entries(str(access),
                                           entrez_attribute,
                                           attribute_suffix,
                                           args.email,
                                           args.number)
                if success != []:
                    sys.stderr.write("\nRetrieval was successful for isolate: " + str(access) + "\n")
                    successful_accessions.add(success)
    sys.exit(0)

if __name__ == '__main__':
    main()
