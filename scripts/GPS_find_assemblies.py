
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script partitions GPS isolates into those with available assemblies, and those without. This is important for determing which pathway they follow in the snakemake pipeline
"""
from Bio import Entrez
from joblib import Parallel, delayed
import pandas as pd
import sys
from tqdm import tqdm

sys.path.insert(0, "..")
from BacQuerya_processing.secrets import ENTREZ_API_KEY

def get_options():

    import argparse

    parser = argparse.ArgumentParser(prog='extract_entrez_reads')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("--metadata-file",
                        dest="metadata_csv",
                        required=True,
                        help='GPS metadata file from which accessions will be extracted',
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

def partition_isolates(cleaned_accession,
                       email,
                       number):
    """Download the biosample accession for accessions of interest using Biopython Entrez"""
    Entrez.email = email
    Entrez.api_key = ENTREZ_API_KEY
    handle = Entrez.read(Entrez.esearch(db="assembly", term=cleaned_accession, retmax = number))
    assembly_ids = handle['IdList']
    read_accessions = []
    try:
        esummary_handle = Entrez.esummary(db="assembly", id=assembly_ids[0], report="full")
        esummary_record = Entrez.read(esummary_handle, validate = False)
        biosample_id = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]["BioSampleAccn"]
        with open("GPS_assembly_biosamples.txt", "a") as outFile:
            outFile.write(biosample_id + "\n")
    except IndexError:
       # there is no assembly available
       with open("GPS_read_laneIDs.txt", "a") as outERS:
            outERS.write(cleaned_accession + "\n")

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()
    GPS_metadata = pd.read_csv(args.metadata_csv)
    sample_list = list(GPS_metadata["Lane_Id"])
    sample_ERS = list(GPS_metadata["ERS"])
    cleaned_accessions = []
    cleaned_ERS = []
    for access in range(len(sample_list)):
        if not sample_list[access] == "":
            cleaned_accessions.append(sample_list[access])
            cleaned_ERS.append(sample_ERS[access])
    del GPS_metadata
    sys.stderr.write("\nPartitioning " + str(len(cleaned_accessions)) + " isolates\n")
    job_list = [
        cleaned_accessions[i:i + args.n_cpu] for i in range(0, len(cleaned_accessions), args.n_cpu)
    ]
    non_assembly_ids = []
    for job in tqdm(job_list):
        Parallel(n_jobs=args.n_cpu)(delayed(partition_isolates)(str(access),
                                                                args.email,
                                                                args.number) for access in job)
    with open("GPS_read_laneIDs.txt", "r") as inERS:
        read_LaneID = inERS.read().splitlines()
    read_ERS = []
    for read in read_LaneID:
        index = sample_list.index(read)
        read_ERS.append(sample_ERS[index])
    sys.stderr.write("\n" + str(len(read_ERS)) + " read accessions, " + str(len(cleaned_accessions) - len(read_ERS)) + " assembly accessions."+"\n")
    with open("GPS_read_ERS.txt", "w") as outFile:
        outFile.write("\n".join(read_ERS))
    sys.exit(0)

if __name__ == '__main__':
    main()
