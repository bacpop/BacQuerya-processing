#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract attributes from isolate-specific assembly stats and constuct a JSON for all attributes in all isolates.
"""
from assembly_stats import read_genome, calculate_stats
from Bio import Entrez
from joblib import Parallel, delayed
import json
import glob
import os
import re
import sys
from tqdm import tqdm
import xmltodict

from BacQuerya_processing.secrets import ENTREZ_API_KEY
from BacQuerya_processing.extract_read_metadata import standardise_species

def get_options():

    import argparse

    description = 'Extract assembly stats'
    parser = argparse.ArgumentParser(description=description,
                                        prog='extract_assembly_stats')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-a",
                        "--assemblies",
                        dest="assemblies",
                        required=True,
                        help='directory of assembly reports',
                        type=str)
    io_opts.add_argument("-g",
                        "--genomes",
                        dest="genomes",
                        required=True,
                        help='directory of genomic sequences',
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
                        help="output file for json of isolates",
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
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    io_opts.add_argument("--GPS",
                        dest="GPS",
                        required=False,
                        help="specify if we are indexing GPS data",
                        action='store_true')
    io_opts.add_argument("--supplementary-metdata",
                        dest="supplementary_metadataJSON",
                        required=False,
                        help="JSON file of supplementary metadata",
                        default=None,
                        type=str)
    args = parser.parse_args()
    return (args)

def get_biosample_metadata(biosample_accession,
                           email):
    "Use Biopython entrez to download metadata found in NCBI BioSample"
    Entrez.email = email
    Entrez.api_key = ENTREZ_API_KEY
    handle = Entrez.read(Entrez.esearch(db="biosample", term=biosample_accession, retmax = 10))
    assembly_ids = handle['IdList']
    try:
        esummary_handle = Entrez.esummary(db="biosample", id=assembly_ids[0], report="full")
        esummary_record = Entrez.read(esummary_handle, validate = False)
    except:
        sys.stderr.write("\nError with:" + biosample_accession + "\n")
        return None
    biosample_identifiers = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]["Identifiers"]
    biosample_metadata_html = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]["SampleData"]
    biosample_metadata_dict = xmltodict.parse(biosample_metadata_html)
    biosample_metadata = {"BioSample_PublicationDate": biosample_metadata_dict["BioSample"]["@publication_date"],
                          "BioSample_LastUpdate": biosample_metadata_dict["BioSample"]["@last_update"],
                          "BioSample_SubmissionDate": biosample_metadata_dict["BioSample"]["@submission_date"],
                          "BioSample_Owner": biosample_metadata_dict["BioSample"]["Owner"]["Name"],
                          "BioSample_Status": biosample_metadata_dict["BioSample"]["Status"]["@status"]}
    try:
        attributes = biosample_metadata_dict["BioSample"]["Attributes"]["Attribute"]
    except TypeError:
        return None
    if isinstance(attributes, list):
        for attr in range(len(attributes)):
            if attributes[attr]["@attribute_name"] == "INSDC center name":
                biosample_metadata.update({"BioSample_INSDCCenterName" : attributes[attr]["#text"]})
            if attributes[attr]["@attribute_name"] == "collection_date":
                biosample_metadata.update({"BioSample_CollectionDate" : attributes[attr]["#text"]})
            if attributes[attr]["@attribute_name"] == "geographic location (country and/or sea)":
                biosample_metadata.update({"BioSample_CollectionLocation" : attributes[attr]["#text"]})
            if attributes[attr]["@attribute_name"] == "host health state":
                biosample_metadata.update({"BioSample_HostHealthState" : attributes[attr]["#text"]})
            if attributes[attr]["@attribute_name"] == "isolation_source":
                biosample_metadata.update({"BioSample_IsolationSource" : attributes[attr]["#text"]})
            if attributes[attr]["@attribute_name"] == "serovar":
                biosample_metadata.update({"BioSample_SeroVar" : attributes[attr]["#text"]})
            if attributes[attr]["@attribute_name"] == "specific_host":
                biosample_metadata.update({"BioSample_SpecificHost" : attributes[attr]["#text"]})
    return biosample_metadata

def calculate_assembly_stats(genomeFile):
    contig_lens, scaffold_lens, gc_cont = read_genome(genomeFile)
    contig_stats = calculate_stats(contig_lens, gc_cont)
    scaffold_stats = calculate_stats(scaffold_lens, gc_cont)
    return contig_stats, scaffold_stats

def assembly_to_JSON(assigned_index,
                     genome_dir,
                     GPS,
                     supplementary_metadataJSON):
    """Use assembly stats to extract information for elasticsearch indexing"""
    index_no = assigned_index['isolate_index']
    assembly_file = assigned_index['assembly file']
    isolate_name = os.path.basename(assembly_file).replace("_assembly_stats.txt", "")
    genome_file = os.path.join(genome_dir, isolate_name + ".fna")
    contig_stats, scaffold_stats = calculate_assembly_stats(genome_file)
    with open(assembly_file, "r") as f:
        assembly_features = f.read().split("\n")
    with open(assembly_file.replace("_assembly_stats.txt", "_additionalAssemblyStats.txt"), "r") as a:
        assemblyURL = json.loads(a.read())
    assembly_dict = {"isolateName": isolate_name.replace("_", " "),
                     "isolateNameUnderscore": isolate_name,
                     "isolate_index": index_no,
                     "contig_stats": contig_stats,
                     "scaffold_stats": scaffold_stats,
                     "sequenceURL": assemblyURL["sequenceURL"],
                     "accession": assemblyURL["accession"]}
    for line in assembly_features:
        try:
            attribute = re.search('# (.*?):', line).group(1).replace(" ", "_")
            feature_dict = {attribute : line.split(":")[-1].strip().replace("_", " ")}
            if "Organism_name" in feature_dict:
                feature_dict["Organism_name"] = standardise_species(feature_dict["Organism_name"])
            assembly_dict.update(feature_dict)
        except AttributeError:
            pass
    if GPS:
        for accession, supplement in supplementary_metadataJSON.items():
            if supplement["Lane_Id"] == assembly_dict["isolateNameUnderscore"]:
                assembly_dict.update(supplement)
    return assembly_dict

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    output_file = os.path.join(args.output_dir, "isolateAssemblyAttributes.json")
    biosampleKeys = os.path.join(args.output_dir, "biosampleIsolatePairs.json")
    isolateKeys = os.path.join(args.output_dir, "indexIsolatePairs.json")
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    assembly_reports = glob.glob(args.assemblies + '/*_assembly_stats.txt')
    indexed_assemblies = []
    with open(args.index_file, "r") as indexFile:
        indexNoDict = json.loads(indexFile.read())
    index_no = int(indexNoDict["isolateIndexNo"])
    # if previous isolate_kv pairs exist use that to prevent duplication of isolates in index
    previousRunFile = os.path.join(args.previous_dir, isolateKeys)
    if os.path.exists(previousRunFile):
        with open(previousRunFile) as prevKeys:
            indexedIsolateDict = json.loads(prevKeys.read())
    else:
        indexedIsolateDict = {}
    for assembly in assembly_reports:
        isol_label = os.path.basename(assembly).replace("_assembly_stats.txt", "")
        if isol_label in indexedIsolateDict.keys():
            index_no = indexedIsolateDict[isol_label]
        else:
            index_no = index_no
        assigned_index = {"isolate_index": index_no,
                          "assembly file": assembly}
        index_no += 1
        indexed_assemblies.append(assigned_index)
        indexedIsolateDict.update({isol_label: index_no})
    # import the GPS metadata JSON if needed
    if args.GPS:
        with open(args.supplementary_metadataJSON, "r") as metaFile:
           supplementary_metadataJSON = json.loads(metaFile.read())
    else:
        supplementary_metadataJSON = None
    sys.stderr.write('\nConverting assembly stat files to JSON\n')
    job_list = [
        indexed_assemblies[i:i + args.n_cpu] for i in range(0, len(indexed_assemblies), args.n_cpu)
    ]
    # parrallelise assembly feature extraction
    all_features = []
    for job in tqdm(job_list):
        features = Parallel(n_jobs=args.n_cpu)(delayed(assembly_to_JSON)(assem,
                                                                         args.genomes,
                                                                         args.GPS,
                                                                         supplementary_metadataJSON) for assem in job)
        all_features += features
    # get label biosample pairs for isolate URL
    sys.stderr.write('\nWriting label : BioSample pairs and extracting BioSample metadata\n')
    labelBiosampleDict = {}
    for feat in tqdm(all_features):
        labelBiosampleDict.update({feat["isolateNameUnderscore"] : feat["BioSample"]})
        biosample_data = get_biosample_metadata(feat["BioSample"], args.email)
        if biosample_data:
            feat.update(biosample_data)
    with open(os.path.join(output_file), "w") as a:
        a.write(json.dumps({"information":all_features}))
    # output isolateName and assigned index k,v pairs
    with open(isolateKeys, "w") as a:
        a.write(json.dumps(indexedIsolateDict))
    # output isolateName and biosample accession k,v pairs
    with open(biosampleKeys, "w") as b:
        b.write(json.dumps(labelBiosampleDict))
    # update isolate index number for subsequent runs
    indexNoDict["isolateIndexNo"] = index_no
    with open(args.index_file, "w") as indexFile:
        indexFile.write(json.dumps(indexNoDict))
    sys.exit(0)

if __name__ == '__main__':
    main()
