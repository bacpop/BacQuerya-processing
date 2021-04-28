import json
import os
import subprocess
import sys

def get_options():

    import argparse

    description = 'Merge the outputa of previous and current snakemake runs'
    parser = argparse.ArgumentParser(description=description,
                                        prog='merge_runs')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("--ncbi-metadata",
                        dest="ncbi_metadata",
                        required=True,
                        help='isolate JSON output by extract_assembly_stats for current run',
                        type=str)
    io_opts.add_argument("--geneMetadataDir",
                        dest="geneMetadataDir",
                        required=True,
                        help="JSON of gene metadata output by extract_genes for current run",
                        type=str)
    io_opts.add_argument("--panarooOutputDir",
                        dest="panarooOutputDir",
                        required=True,
                        help="Directory output by panaroo for current run",
                        type=str)
    io_opts.add_argument("--accessionFile",
                        dest="accessionFile",
                        required=True,
                        help="File of accession IDs used in the current run",
                        type=str)
    io_opts.add_argument("--previous-run",
                        dest="prev_run",
                        required=True,
                        help="directory of outputs containing information from previous runs",
                        type=str)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    args = parser.parse_args()
    return (args)

def mergeNCBIMetadata(current_metadata_file, prev_dir):
    with open(current_metadata_file, "r") as currentFile:
        currentIsolateJSON = currentFile.read()
    previous_metadata_file = os.path.join(prev_dir, os.path.basename(current_metadata_file))
    with open(previous_metadata_file, "r") as previousFile:
        previousIsolateJSON = previousFile.read()
    currentIsolateDicts = json.loads(currentIsolateJSON)["information"]
    updatedIsolateDict = json.loads(previousIsolateJSON)
    updatedIsolateDict["information"] += currentIsolateDicts
    with open(previous_metadata_file, "w") as updatedFile:
        updatedFile.write(json.dumps(updatedIsolateDict))

def mergeNCBIKVPairs(current_KVPairs, current_biosamplePairs, prev_dir):
    # update isolate key value json
    with open(current_KVPairs, "r") as currentKVFile:
        currentKVJSON = currentKVFile.read()
    previous_KVPairs = os.path.join(prev_dir, os.path.basename(current_KVPairs))
    with open(previous_KVPairs, "r") as previousKVFile:
        previousKVJSON = previousKVFile.read()
    currentKVDict = json.loads(currentKVJSON)
    current_isolateIndex = max(currentKVDict.values())
    updatedKVDict = json.loads(previousKVJSON)
    updatedKVDict.update(currentKVDict)
    with open(previous_KVPairs, "w") as updatedKVFile:
        updatedKVFile.write(json.dumps(updatedKVDict))
    # update biosample isolate pairs json
    with open(current_biosamplePairs, "r") as currentBiosampleFile:
        currentBiosampleJSON = currentBiosampleFile.read()
    previous_biosamplePairs = os.path.join(prev_dir, os.path.basename(current_biosamplePairs))
    with open(previous_biosamplePairs, "r") as previousBiosampleFile:
        previousBiosampleJSON = previousBiosampleFile.read()
    currentBiosampltDict = json.loads(currentBiosampleJSON)
    updatedBiosampleDict = json.loads(previousBiosampleJSON)
    updatedBiosampleDict.update(currentBiosampltDict)
    with open(previous_biosamplePairs, "w") as updatedBiosampleFile:
        updatedBiosampleFile.write(json.dumps(updatedBiosampleDict))
    return current_isolateIndex

def mergeGeneMetadata(current_geneDir, prev_dir):
    # update gene metadata json
    current_metadataFile = os.path.join(current_geneDir, "annotatedNodes.json")
    with open(current_metadataFile, "r") as currentGeneFile:
        currentGeneJSON = currentGeneFile.read()
    previous_metadataFile = os.path.join(prev_dir, current_metadataFile)
    with open(previous_metadataFile, "r") as previousGeneFile:
        previousGeneJSON = previousGeneFile.read()
    currentGeneDicts = json.loads(currentGeneJSON)["information"]
    updatedGeneDict = json.loads(previousGeneJSON)
    updatedGeneDict["information"] += currentGeneDicts
    with open(previous_metadataFile, "w") as updatedGeneFile:
       updatedGeneFile.write(json.dumps(updatedGeneDict))
    # update gene: index key value pairs
    current_panarooPairs = os.path.join(current_geneDir, "panarooPairs.json")
    with open(current_panarooPairs, "r") as currentKVFile:
        currentKVJSON = currentKVFile.read()
    previous_panarooPairs = os.path.join(prev_dir, current_panarooPairs)
    with open(previous_panarooPairs, "r") as previousKVFile:
        previousKVJSON = previousKVFile.read()
    currentKVDict = json.loads(currentKVJSON)
    currentGeneIndex = max(currentKVDict.values())
    updatedKVDict = json.loads(previousKVJSON)
    updatedKVDict.update(currentKVDict)
    with open(previous_panarooPairs, "w") as updatedKVFile:
       updatedKVFile.write(json.dumps(updatedKVDict))
    return currentGeneIndex

def mergeAccessionIDs(current_accessionFile, prev_dir):
    with open(current_accessionFile, "r") as currentFile:
        current_accessionList = currentFile.read().splitlines()
    previous_accessionFile = os.path.join(prev_dir, os.path.basename(current_accessionFile))
    with open(previous_accessionFile, "r") as previousFile:
        previous_accessionList = previousFile.read().splitlines()
    updated_accessionSet = set(previous_accessionList)
    for access in current_accessionList:
        updated_accessionSet.add(access)
    updated_accessionList = list(updated_accessionSet)
    with open(previous_accessionFile, "w") as updatedFile:
        updatedFile.write("\n".join(updated_accessionList))

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(args.prev_run):
        sys.stderr.write("\nCopying current run data into " + args.prev_run + "\n")
        subprocess_command = "mkdir "
        subprocess_command += args.prev_run
        subprocess_command += " && cp -r "
        subprocess_command += args.ncbi_metadata + " "
        subprocess_command += args.geneMetadataDir + " "
        subprocess_command += args.panarooOutputDir + " "
        subprocess_command += args.accessionFile + " "
        subprocess_command += " previous_run"
        subprocess.run(subprocess_command, shell=True, check=True)
    else:
        # merge metadata for isolates found in NCBI for current and previous runs
        sys.stderr.write("\nMerging current and previous NCBI isolate metadata\n")
        currentIolateMetadata = os.path.join(args.ncbi_metadata, "isolateAssemblyAttributes.json")
        mergeNCBIMetadata(currentIolateMetadata,
                          args.prev_run)
        sys.stderr.write("\nDone\n")
        # merge key value pairs for current and previous runs
        sys.stderr.write("\nMerging current and previous NCBI key value pairs\n")
        current_indexIsolatePairs = os.path.join(args.ncbi_metadata, "indexIsolatePairs.json")
        current_biosampleIsolatePairs = os.path.join(args.ncbi_metadata, "biosampleIsolatePairs.json")
        currentIsolateIndex = mergeNCBIKVPairs(current_indexIsolatePairs,
                                               current_biosampleIsolatePairs,
                                               args.prev_run)
        sys.stderr.write("\nDone\n")
        # merge gene metadata for current and previous runs
        sys.stderr.write("\nMerging current and previous gene metadata\n")
        currentGeneIndex = mergeGeneMetadata(args.geneMetadataDir,
                                             args.prev_run)
        # merge panaroo output for current and previous runs
        sys.stderr.write("\nMerging current and previous panaroo outputs\n")
        previous_panarooOutput = os.path.join(args.prev_run, args.panarooOutputDir)
        panarooCommand = "panaroo-merge -d " + args.panarooOutputDir + " " + previous_panarooOutput + " -o " + previous_panarooOutput + " -t " + str(args.n_cpu)
        subprocess.run(panarooCommand, shell=True, check=True)
        sys.stderr.write("\nDone\n")
        # merge accession IDs for current and previous runs
        sys.stderr.write("\nMerging current and previous accession IDs\n")
        mergeAccessionIDs(args.accessionFile,
                          args.prev_run)
        sys.stderr.write("\nDone\n")
    sys.exit(0)
