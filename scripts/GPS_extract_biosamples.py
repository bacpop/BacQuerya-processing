
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses Biopython ENTREZ to retrieve and download information of interest available through NCBI.
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

    description = 'Download reads of interest by accession ID'
    parser = argparse.ArgumentParser(description=description,
                                        prog='extract_entrez_reads')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("--metadata-file",
                        dest="metadata_csv",
                        required=True,
                        help='GPS metadata file from which accessions will be extracted',
                        type=str)
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_file",
                        required=True,
                        help="name of file listing GPS biosample accessions",
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
                          number):
    """Download the biosample accession for accessions of interest using Biopython Entrez"""
    successful_accessions = []
    Entrez.email = email
    Entrez.api_key = ENTREZ_API_KEY
    handle = Entrez.esearch(db="sra", term=cleaned_accession, retmax = number)
    result = Entrez.read(handle)
    try:
        biosample_id = result['IdList'][0]
        print(biosample_id)
        successful_accessions.append(biosample_id)
    except:
        # issue with the request. Re-requesting often solves
        sys.stderr.write("\nRequest failed, the following accession will be re-requested: " + cleaned_accession + "\n")
    return successful_accessions

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()
    GPS_metadata = pd.read_csv(args.metadata_csv)
    sample_list = list(GPS_metadata["ERS"])
    cleaned_accessions = []
    for access in range(len(sample_list)):
        if not sample_list[access] == "":
            cleaned_accessions.append(sample_list[access])
    del GPS_metadata
    # set wd to output_dir for urllib
    sys.stderr.write("\nDownloading biosamples for " + str(len(cleaned_accessions)) + " isolates\n")
    job_list = [
        cleaned_accessions[i:i + args.n_cpu] for i in range(0, len(cleaned_accessions), args.n_cpu)
    ]
    successful_accessions = []
    for job in tqdm(job_list):
        successful_accessions += Parallel(n_jobs=args.n_cpu)(delayed(download_SRA_metadata)(access,
                                                                                            args.email,
                                                                                            args.number) for access in job)
    # ensure available attributes are downloaded for all accessions of interest
    successful_accessions = [success for row in successful_accessions for success in row]
    successful_accessions = set(successful_accessions)
    while len(successful_accessions) != len(cleaned_accessions):
        for access in cleaned_accessions:
            if not access in successful_accessions:
                sys.stderr.write("\nRerequesting isolate: " + str(access) + "\n")
                success = download_SRA_metadata(str(access),
                                                args.email,
                                                args.number)
                if success != []:
                    sys.stderr.write("\nRetrieval was successful for isolate: " + str(access) + "\n")
                    successful_accessions.add(success[0])
    sys.exit(0)

if __name__ == '__main__':
    main()
