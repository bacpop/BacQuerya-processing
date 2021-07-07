#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses Biopython ENTREZ to retrieve and download metadata for reads of interest and output as a json.
"""
from Bio import Entrez
from joblib import Parallel, delayed
import json
import os
import requests
from shutil import copyfileobj, rmtree
import ssl
import subprocess
import sys
import tempfile
from tqdm import tqdm
import xmltodict

from BacQuerya_processing.extract_assembly_stats import get_biosample_metadata, calculate_assembly_stats
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
    io_opts.add_argument("-m",
                        "--isolate-metadata",
                        dest="isolate_metadata",
                        required=True,
                        help='Directory output by extract_assembly_stats.py',
                        type=str)
    io_opts.add_argument("-i",
                        "--index-file",
                        dest="index_file",
                        required=True,
                        help="JSON file containing integer value to start index from",
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
    io_opts.add_argument("--previous-run",
                        dest="previous_dir",
                        required=True,
                        help="directory name of previous snakemake outputs",
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
    io_opts.add_argument("--GPS",
                        dest="GPS",
                        required=False,
                        help="specify if we are indexing GPS data",
                        action='store_true')
    io_opts.add_argument("--ESC",
                        dest="ESC",
                        required=False,
                        help="specify if we are indexing ESC data",
                        action='store_true')
    io_opts.add_argument("--supplementary-metdata",
                        dest="supplementary_metadataJSON",
                        required=False,
                        help="JSON file of supplementary metadata",
                        default=None,
                        type=str)
    io_opts.add_argument("--assembly-url",
                        dest="assemblyURLs",
                        required=False,
                        help="JSON file of biosample assembly pairs output by scripts/GPS_661kassembly_links.py",
                        default=None,
                        type=str)
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

def get_ENA_accession(Id):
    """Search for an Id in the ENA if a BioSample, run or read accession is not provided"""
    apiURL = "https://www.ebi.ac.uk/ebisearch/ws/rest/sra-sample?query={}&fields=acc,description,name&size=100".format(Id.replace("#", "%23"))
    urlResponse = requests.get(apiURL)
    accession_metadata = dict(xmltodict.parse(urlResponse.text))
    try:
        if accession_metadata["result"]["hitCount"] == "0":
            apiURL = "https://www.ebi.ac.uk/ebisearch/ws/rest/wgs_masters?query={}&fields=acc,description,name&size=100".format(Id.replace("#", "%23"))
            urlResponse = requests.get(apiURL)
            accession_metadata = dict(xmltodict.parse(urlResponse.text))
            if accession_metadata["result"]["hitCount"] == "0":
                return None
            else:
                try:
                    access = accession_metadata["result"]["entries"]['entry']["@id"]
                except:
                    return None
                txtURL = "https://www.ebi.ac.uk/ena/browser/api/embl/{}?lineLimit=1000".format(access)
                urlResponse = requests.get(txtURL).text.splitlines()
                for resp in urlResponse:
                    if "BioSample" in resp:
                        biosample = resp.split(" ")[4].replace(".", "")
                        return biosample
        else:
            access = accession_metadata["result"]["entries"]['entry']["@id"]
            return access
    except:
        return None

def standardise_species(organism_name):
    """Standardise the species name for the BacQuerya website"""
    species_dict = {"Streptococcus pneumoniae": "Streptococcus pneumoniae",
                "S. pneumoniae": "Streptococcus pneumoniae",
                "S.pneumoniae": "Streptococcus pneumoniae",
                "Escherichia coli": "Escherichia coli",
                "E. coli": "Escherichia coli",
                "E.coli": "Escherichia coli",
                "Staphylococcus aureus": "Staphylococcus aureus",
                "S. aureus": "Staphylococcus aureus",
                "S.aureus": "Staphylococcus aureus"}
    for key, value in species_dict.items():
        if value in organism_name:
            return key

def download_ENA_metadata(accession_dict,
                          temp_dir,
                          GPS,
                          ESC,
                          supplementary_metadataJSON,
                          assemblyURLs):
    cleaned_accession = accession_dict["input_accession"]
    queried_accession = cleaned_accession
    index_no = accession_dict["isolate_index"]
    apiURL = "https://www.ebi.ac.uk/ena/browser/api/xml/" + cleaned_accession
    #try:
    urlResponse = requests.get(apiURL)
    accession_metadata = dict(xmltodict.parse(urlResponse.text))
    if "ErrorDetails" in accession_metadata.keys():
        identifier = get_ENA_accession(cleaned_accession)
        if identifier:
            apiURL = "https://www.ebi.ac.uk/ena/browser/api/xml/" + identifier
            urlResponse = requests.get(apiURL)
            accession_metadata = dict(xmltodict.parse(urlResponse.text))
        else:
            with open("ENA_READS_ERROR.txt", "a") as suppressed:
                suppressed.write(cleaned_accession + "\n")
            return None
    if "RUN_SET" in accession_metadata.keys():
        # we're looking at the run accession,, need to extract the read accession
        temp_metadata = accession_metadata["RUN_SET"]["RUN"]
        run_accession = cleaned_accession
        if isinstance(temp_metadata["RUN_LINKS"]["RUN_LINK"], list):
            for elem in temp_metadata["RUN_LINKS"]["RUN_LINK"]:
                if elem["XREF_LINK"]["DB"] == "ENA-SAMPLE":
                    read_accession = elem["XREF_LINK"]["ID"]
        else:
            if "XREF_LINK" in temp_metadata["RUN_LINKS"]["RUN_LINK"] and temp_metadata["RUN_LINKS"]["RUN_LINK"]["XREF_LINK"]["DB"] == "ENA-SAMPLE":
                read_accession = temp_metadata["RUN_LINKS"]["RUN_LINK"]["XREF_LINK"]["ID"]
        apiURL = "https://www.ebi.ac.uk/ena/browser/api/xml/" + read_accession
        urlResponse = requests.get(apiURL)
        accession_metadata = dict(xmltodict.parse(urlResponse.text))
    elif "SAMPLE_SET" in accession_metadata.keys():
        cleaned_accession = cleaned_accession
        run_accession = ""
    try:
        accession_metadata = accession_metadata["SAMPLE_SET"]
    except:
        with open("ENA_READS_ERROR.txt", "a") as error:
            error.write(queried_accession + "\n")
        return None
    # extract biosample id for isolate
    try:
        biosample_id = accession_metadata["SAMPLE"]["IDENTIFIERS"]["EXTERNAL_ID"]["#text"]
    except TypeError:
        for ext in accession_metadata["SAMPLE"]["IDENTIFIERS"]["EXTERNAL_ID"]:
            if ext["@namespace"] == "BioSample":
                biosample_id = ext["#text"]
    assembly_stats = False
    for link in accession_metadata["SAMPLE"]["SAMPLE_LINKS"]["SAMPLE_LINK"]:
        if link["XREF_LINK"]["DB"] == "ENA-RUN":
            run_accession = link["XREF_LINK"]["ID"]
        if link["XREF_LINK"]["DB"] == "ENA-FASTQ-FILES":
            fastqTable = requests.get(link["XREF_LINK"]["ID"]).text.split("\t")
            try:
                fastqLinks = fastqTable[4].split(";")
                for fql in range(len(fastqLinks)):
                    fastqLinks[fql] = "https://" + fastqLinks[fql]
            except:
                fastqLinks = []
            # if an assembly is available in the ENA, add the link to the sequence links
            if assemblyURLs:
                if biosample_id in assemblyURLs:
                    # if a Blackwell assembly is available, we need to retrieve it and calculate assembly statistics
                    genome_representation = "full"
                    assemblyLink = assemblyURLs[biosample_id]
                    # download the assembly file and save to a temp directory
                    #ssl._create_default_https_context = ssl._create_unverified_context
                    assemblyFile = os.path.join(temp_dir, os.path.basename(assemblyLink))
                    #with urlopen(assemblyLink) as in_stream, open(assemblyFile, 'wb') as out_file:
                        #copyfileobj(in_stream, out_file)
                    subprocess.run("wget -q -O " + assemblyFile + " " + assemblyLink.replace("http", "ftp"), shell = True, check = True)
                    #urllib.request.urlretrieve(assemblyLink.replace("http", "ftp"), assemblyFile)
                    # unzip the assembly file if necessary
                    try:
                        if ".gz" in assemblyFile:
                            subprocess.run("gunzip " + assemblyFile, shell=True, check=True)
                        # calculate assembly statistics
                        contig_stats, scaffold_stats = calculate_assembly_stats(assemblyFile.replace(".gz", ""))
                        assembly_stats = True
                    except:
                        assembly_stats = False
                    fastqLinks.append(assemblyLink)
                else:
                    genome_representation = "reads"
            accession_metadata.update({"ENA-FASTQ-FILES" : fastqLinks})
    for attribute in accession_metadata["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]:
        if attribute["TAG"] == "ENA-FIRST-PUBLIC":
            submission_date = attribute["VALUE"]
    try:
        submitter = accession_metadata["SAMPLE"]["IDENTIFIERS"]["SUBMITTER_ID"]["@namespace"]
    except:
        try:
            submitter = accession_metadata["SAMPLE"]["IDENTIFIERS"]["SUBMITTER_ID"]["namespace"]
        except KeyError:
            submitter = None
    # if we are indexing the GPS data, we need to set the isolate name as the lane_id
    isolateName = cleaned_accession
    if GPS:
        try:
            GPS_metadata = supplementary_metadataJSON[run_accession]
            isolateName = GPS_metadata["Lane_Id"].replace("_", " ")
        except KeyError:
            for accession, supplement in supplementary_metadataJSON.items():
                if "ERS" in supplement.keys() and supplement["ERS"] == cleaned_accession:
                    GPS_metadata = supplement
                    isolateName = GPS_metadata["Lane_Id"].replace("_", " ")
    if ESC:
        ESC_metadata = supplementary_metadataJSON[queried_accession]
        isolateName = cleaned_accession
    metadata = {"isolateName" : isolateName,
                "accession" : cleaned_accession,
                "read_accession" : cleaned_accession,
                "isolate_index" : index_no,
                "Genome_representation" : genome_representation,
                "SubmissionDate" : submission_date,
                "Organism_name" : standardise_species(accession_metadata["SAMPLE"]["SAMPLE_NAME"]["SCIENTIFIC_NAME"]),
                "Taxid" : accession_metadata["SAMPLE"]["SAMPLE_NAME"]["TAXON_ID"],
                "BioSample" : biosample_id,
                "source" : "ENA",
                "sequenceURL" : accession_metadata["ENA-FASTQ-FILES"],
                "allAttributes" : json.dumps(accession_metadata["SAMPLE"])}
    if GPS:
        metadata.update(GPS_metadata)
    if ESC:
        if not "Year" in ESC_metadata and "BioSample_CollectionDate" in ESC_metadata:
            ESC_metadata["Year"] = int(ESC_metadata["BioSample_CollectionDate"])
            print(ESC_metadata["Year"])
        metadata.update(ESC_metadata)
    if submitter:
        metadata.update({"Submitter": submitter})
    # add assembly stats to isolate metadata if it is defined
    if assembly_stats:
        metadata["contig_stats"] = contig_stats
        metadata["scaffold_stats"] = scaffold_stats
    # output isolate: index pairs
    indexIsolatePair = {isolateName: index_no}
    if not run_accession == "":
        metadata.update({"run_accession" : run_accession})
    return [fastqLinks, metadata, indexIsolatePair]
    #except:
       # sys.stderr.write("Request failed for the following accession: " + cleaned_accession)
       # with open("ENA_READS_ERROR.txt", "a") as failed:
         #   failed.write(cleaned_accession + "\n")
        #return None

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
    with open(args.index_file, "r") as indexFile:
        indexNoDict = json.loads(indexFile.read())
    index_no = int(indexNoDict["isolateIndexNo"])
    # if previous isolate_kv pairs exist use that to prevent duplication of isolates in index
    previousRunFile = os.path.join(args.previous_dir, args.output_dir, "indexIsolatePairs.json")
    if os.path.exists(previousRunFile):
        with open(previousRunFile) as prevKeys:
            indexIsolateDict = json.loads(prevKeys.read())
    else:
        indexIsolateDict = {}
    for access in range(len(sample_list)):
        if not sample_list[access] == "":
            if sample_list[access] in indexIsolateDict.keys():
                index_no = indexIsolateDict[sample_list[access]]
            else:
                index_no = index_no
            assigned_index = {"isolate_index": index_no,
                              "input_accession": sample_list[access]}
            index_no += 1
            indexed_accessions.append(assigned_index)
    job_list = [
        indexed_accessions[i:i + args.n_cpu] for i in range(0, len(indexed_accessions), args.n_cpu)
    ]
    # import the GPS metadata JSON if needed
    if args.GPS or args.ESC:
        with open(args.supplementary_metadataJSON, "r") as metaFile:
           supplementary_metadataJSON = json.loads(metaFile.read())
    else:
        supplementary_metadataJSON = None
    # import file of biosample assembly links k, v pairs
    if args.assemblyURLs:
        with open(args.assemblyURLs, "r") as linksFile:
           assemblyURLs = json.loads(linksFile.read())
    else:
        assemblyURLs = None
    if args.read_source == "ena":
        # need a tempdir to dowload the Blackwell assemblies
        temp_dir = tempfile.mkdtemp(dir=args.output_dir)
        access_data = []
        for job in tqdm(job_list):
            access_data +=  Parallel(n_jobs=args.n_cpu)(delayed(download_ENA_metadata)(access,
                                                                                       temp_dir,
                                                                                       args.GPS,
                                                                                       args.ESC,
                                                                                       supplementary_metadataJSON,
                                                                                       assemblyURLs) for access in job)
        #access_data = [link for row in access_data for link in row]
        fastq_links = []
        metadata = []
        for line in access_data:
            if line:
                fastq_links += line[0]
                metadata.append(line[1])
                indexIsolateDict.update(line[2])
        # get BioSample Metadata
        sys.stderr.write('\nDownloading BioSample metadata\n')
        #### need to get assembly accessions too if they are present for the GPS data
        # also need to store biosample accessions for read data
        biosample_pairs = {}
        for metadata_line in tqdm(metadata):
            biosample_metadata = get_biosample_metadata(metadata_line["BioSample"], args.email)
            if biosample_metadata:
                metadata_line.update(biosample_metadata)
                biosample_pairs[metadata_line["isolateName"]] = metadata_line["BioSample"]
        # write out list of run accessions
        with open(os.path.join(args.output_dir, "fastq_links.txt"), "w") as r:
            r.write("\n".join(fastq_links))
        metadata_json = {"information" : metadata}
        metadata_json = json.dumps(metadata_json).replace("@", "")
        with open(os.path.join(args.output_dir, "isolateReadAttributes.json"), "w") as f:
            f.write(metadata_json)
        with open(os.path.join(args.output_dir, "indexIsolatePairs.json"), "w") as indexPairs:
            indexPairs.write(json.dumps(indexIsolateDict))
        # update isolate index number for subsequent runs
        indexNoDict["isolateIndexNo"] = index_no
        with open(args.index_file, "w") as indexFile:
            indexFile.write(json.dumps(indexNoDict))
        # update isolate biosample label dict
        with open(os.path.join(args.isolate_metadata, "biosampleIsolatePairs.json"), "r") as inBiosample:
            biosample_dict = json.loads(inBiosample.read())
        biosample_dict.update(biosample_pairs)
        with open(os.path.join(args.isolate_metadata, "biosampleIsolatePairs.json"), "w") as outBiosample:
            outBiosample.write(json.dumps(biosample_dict))
        # remove temp dir
        rmtree(temp_dir)
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
