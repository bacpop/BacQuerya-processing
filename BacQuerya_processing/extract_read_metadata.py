#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses Biopython ENTREZ to retrieve and download metadata for reads of interest and output as a json.
"""
from Bio import Entrez
from joblib import Parallel, delayed
import json
import os
import sys
import subprocess
from tqdm import tqdm
import requests
import xmltodict

from BacQuerya_processing.extract_assembly_stats import get_biosample_metadata
from BacQuerya_processing.secrets import ENTREZ_API_KEY

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
    io_opts.add_argument("-r",
                        "--read-source",
                        dest="read_source",
                        required=True,
                        help='Source for read data (SRA or ENA)',
                        choices=['sra', 'ena'],
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


def download_SRA_metadata(cleaned_accession,
                          email,
                          number,
                          output_dir):
    """Download the attribute for accessions of interest using Biopython Entrez"""
    failed_accessions = []
    Entrez.email = email
    Entrez.api_key = ENTREZ_API_KEY
    handle = Entrez.read(Entrez.esearch(db="sra", term=cleaned_accession, retmax = number))
    assembly_ids = handle['IdList']
    try:
        esummary_handle = Entrez.esummary(db="biosample", id=assembly_ids[0], report="full")
        esummary_record = Entrez.read(esummary_handle, validate = False)
        accession_metadata = dict(esummary_record["DocumentSummarySet"]["DocumentSummary"][0])
        accession_metadata["SampleData"] = xmltodict.parse(accession_metadata["SampleData"])
        metadata = {cleaned_accession : accession_metadata}
        metadata_json = json.dumps(metadata).replace("@", "")
        with open(os.path.join(output_dir, cleaned_accession + ".json"), "w") as f:
            f.write(metadata_json)
    except ValueError:
        # the requested attribute is not present
        sys.stderr.write("Metadata is not available for: " + cleaned_accession)
        pass
    except:
        # issue with the request. Re-requesting often solves
        failed_accessions.append(cleaned_accession)
        sys.stderr.write("Request failed, the following accession will be re-requested: " + cleaned_accession)
    return failed_accessions

def download_ENA_metadata(accession_dict,
                          output_dir):
    cleaned_accession = accession_dict["read_accession"]
    index_no = accession_dict["isolate_index"]
    apiURL = "https://www.ebi.ac.uk/ena/browser/api/xml/" + cleaned_accession
    urlResponse = requests.get(apiURL)
    accession_metadata = dict(xmltodict.parse(urlResponse.text))["SAMPLE_SET"]
    try:
        run_accession = ""
        for link in accession_metadata["SAMPLE"]["SAMPLE_LINKS"]["SAMPLE_LINK"]:
            if link["XREF_LINK"]["DB"] == "ENA-RUN":
                run_accession = link["XREF_LINK"]["ID"]
            if link["XREF_LINK"]["DB"] == "ENA-FASTQ-FILES":
                fastqTable = requests.get(link["XREF_LINK"]["ID"]).text.split("\t")
                fastqLinks = fastqTable[4].split(";")
                for fql in range(len(fastqLinks)):
                    fastqLinks[fql] = "https://" + fastqLinks[fql]
                accession_metadata.update({"ENA-FASTQ-FILES" : fastqLinks})
        for attribute in accession_metadata["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]:
            if attribute["TAG"] == "ENA-FIRST-PUBLIC":
                submission_date = attribute["VALUE"]
        try:
            submitter = accession_metadata["SAMPLE"]["IDENTIFIERS"]["SUBMITTER_ID"]["@namespace"]
        except:
            submitter = accession_metadata["SAMPLE"]["IDENTIFIERS"]["SUBMITTER_ID"]["namespace"]
        metadata = {"isolateName" : cleaned_accession,
                    "accession" : cleaned_accession,
                    "read_accession" : cleaned_accession,
                    "isolate_index" : index_no,
                    "Submitter" : submitter,
                    "Genome_representation" : "reads",
                    "Date" : submission_date,
                    "Organism_name" : accession_metadata["SAMPLE"]["SAMPLE_NAME"]["SCIENTIFIC_NAME"],
                    "Taxid" : accession_metadata["SAMPLE"]["SAMPLE_NAME"]["TAXON_ID"],
                    "BioSample" : accession_metadata["SAMPLE"]["IDENTIFIERS"]["EXTERNAL_ID"]["#text"],
                    "source" : "ENA",
                    "sequenceURL" : accession_metadata["ENA-FASTQ-FILES"],
                    "allAttributes" : json.dumps(accession_metadata["SAMPLE"])}
        # output isolate: index pairs
        indexIsolatePair = {cleaned_accession: index_no}
        if not run_accession == "":
            metadata.update({"run_accession" : run_accession})
        metadata_json = json.dumps(metadata).replace("@", "")
    except:
       sys.stderr.write("Request failed for the following accession: " + cleaned_accession)
    return [fastqLinks, metadata, indexIsolatePair]

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    with open(args.accessions, "r") as s:
        sample_accessions = s.read()
    sample_list = sample_accessions.split("\n")
    # assign isolate index number to read metadata
    indexed_accessions = []
    index_no = args.index_no
    for access in range(len(sample_list)):
        if not sample_list[access] == "":
            assigned_index = {"isolate_index": index_no,
                              "read_accession": sample_list[access]}
            index_no += 1
            indexed_accessions.append(assigned_index)
    job_list = [
        indexed_accessions[i:i + args.n_cpu] for i in range(0, len(indexed_accessions), args.n_cpu)
    ]
    indexIsolateDict = {}
    if args.read_source == "ena":
        access_data = []
        for job in tqdm(job_list):
            access_data +=  Parallel(n_jobs=args.n_cpu)(delayed(download_ENA_metadata)(access,
                                                                                       args.output_dir) for access in job)
        #access_data = [link for row in access_data for link in row]
        fastq_links = []
        metadata = []
        for line in access_data:
            fastq_links += line[0]
            metadata.append(line[1])
            indexIsolateDict.update(line[2])
        # get BioSample Metadata
        sys.stderr.write('\nDownload BioSample metadata\n')
        for metadata_line in tqdm(metadata):
            biosample_metadata = get_biosample_metadata(metadata_line["BioSample"], args.email)
            metadata_line.update(biosample_metadata)
        # write out list of run accessions
        with open(os.path.join(args.output_dir, "fastq_links.txt"), "w") as r:
            r.write("\n".join(fastq_links))
        metadata_json = {"information" : metadata}
        metadata_json = json.dumps(metadata_json).replace("@", "")
        with open(os.path.join(args.output_dir, "isolateReadAttributes.json"), "w") as f:
            f.write(metadata_json)
        with open(os.path.join(args.output_dir, "indexIsolatePairs.json"), "w") as indexPairs:
            indexPairs.write(json.dumps(indexIsolateDict))
    if args.read_source == "sra":
        failed_accessions = []
        for job in tqdm(job_list):
            failed_accessions += Parallel(n_jobs=args.n_cpu)(delayed(download_SRA_metadata)(access,
                                                                                            args.email,
                                                                                            args.number,
                                                                                            args.output_dir) for access in job)
        failed_accessions = [failed for row in failed_accessions for failed in row]
        # ensure available attributes are downloaded for all accessions of interest
        for failed_access in failed_accessions:
            failed = download_SRA_metadata(str(failed_access),
                                           args.email,
                                           args.number,
                                           args.output_dir)
            if not len(failed) == 0:
                failed_accessions.append(failed)
    sys.exit(0)

if __name__ == '__main__':
    main()
