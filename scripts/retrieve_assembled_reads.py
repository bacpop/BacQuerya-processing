from joblib import Parallel, delayed
import os
import ssl
import sys
from urllib.request import urlopen
from shutil import copyfileobj
from tqdm import tqdm

def download_read(accession, output_dir):
    if "contigs" in accession:
        try:
            ## currently only downloading assemblies and not read sets for efficiency
            ssl._create_default_https_context = ssl._create_unverified_context
            with urlopen(accession) as in_stream, open(os.path.join(output_dir, os.path.basename(accession)), 'wb') as out_file:
                copyfileobj(in_stream, out_file)
            return [accession]
        except:
            sys.stderr.write("\nRequest failed with: " + accession + "\n")

output_dir = "retrieved_ena_reads"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
with open("fastq_links.txt", "r") as f:
    run_accessions = f.read().splitlines()
cleaned_accessions = []
for line in run_accessions:
    if "contigs" in line:
        cleaned_accessions.append(line)
job_list = [
    cleaned_accessions[i:i + 8] for i in range(0, len(cleaned_accessions), 8)
]
results = []
for job in tqdm(job_list):
        results += Parallel(n_jobs=8)(delayed(download_read)(access,
                                                            output_dir) for access in job)
results = set(results)

while len(results) != len(cleaned_accessions):
    for access in cleaned_accessions:
        if not access in results:
            sys.stderr.write("\nRerequesting: " + str(access) + "\n")
            success = download_read(access,
                                    output_dir)
            if success != []:
                sys.stderr.write("\nRetrieval was successful for: " + str(access) + "\n")
                results.add(success[0])