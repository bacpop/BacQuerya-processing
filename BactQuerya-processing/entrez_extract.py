#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses Biopython ENTREZ to retrieve and download information of interest available through NCBI.
- multiprocess download
"""
from Bio import Entrez
import urllib.request
from tqdm import tqdm
from joblib import Parallel, delayed
import sys
import os

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
                        help='File of sequences to download. One line per entry. Specify accessions of interest or "genus species" for all accessions.',
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
    """Download the genomic attribute for accessions of interest using Biopython Entrez"""
    failed_accessions = []
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
        attribute_link = os.path.join(url,label + entrez_attribute)
        urllib.request.urlretrieve(attribute_link, label + attribute_suffix)
    except ValueError: 
        # the requested attribute is not present
        sys.stderr.write("The requested attribute is not present for: " + cleaned_accession)
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

    entrez_attribute, attribute_suffix = translate_attr(args.attrs)

    # set wd to output_dir for urllib
    os.chdir(args.output_dir)
    job_list = [
        cleaned_accessions[i:i + args.n_cpu] for i in range(0, len(cleaned_accessions), args.n_cpu)
    ]
    for job in tqdm(job_list):
        failed_accessions = Parallel(n_jobs=args.n_cpu)(delayed(download_entries)(access, 
                                                                                entrez_attribute, 
                                                                                attribute_suffix, 
                                                                                args.email, 
                                                                                args.number) for access in job)
    failed_accessions = [failed for row in failed_accessions for failed in row]
    # ensure available attributes are downloaded for all accessions of interest
    for failed_access in failed_accessions:
        failed = download_entries(failed_access,  
                                entrez_attribute, 
                                attribute_suffix,  
                                args.email,
                                args.number)
        if not len(failed) == 0:
            failed_accessions.append(failed)
            print(failed_accessions)
    sys.exit(0)

if __name__ == '__main__':
    main()