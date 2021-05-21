import glob
from joblib import Parallel, delayed
import json
import os
import subprocess
import sys
from tqdm import tqdm

def get_options():

    import argparse

    description = 'Merge the outputa of previous and current snakemake runs'
    parser = argparse.ArgumentParser(description=description,
                                        prog='merge_runs')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("--ncbi-metadata",
                        dest="ncbi_metadata",
                        required=True,
                        help='directory output by extract_assembly_stats for current run',
                        type=str)
    io_opts.add_argument("--geneMetadataDir",
                        dest="geneMetadataDir",
                        required=True,
                        help="directory of gene metadata output by extract_genes for current run",
                        type=str)
    io_opts.add_argument("--alignment-dir",
                        dest="alignment_dir",
                        required=True,
                        help="directory of alignments output by generate_alignments for current run",
                        type=str)
    io_opts.add_argument("--graph-dir",
                        dest="graph_dir",
                        required=True,
                        help="directory of merged panaroo runs",
                        type=str)
    io_opts.add_argument("--accessionFile",
                        dest="accessionFile",
                        required=True,
                        help="file of accession IDs used in the current run",
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
    previous_metadata_file = os.path.join(prev_dir, current_metadata_file)
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
    previous_KVPairs = os.path.join(prev_dir, current_KVPairs)
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
    previous_biosamplePairs = os.path.join(prev_dir, current_biosamplePairs)
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
    updatedGeneDict["information"].update(currentGeneDicts)
    with open(previous_metadataFile, "w") as updatedGeneFile:
       updatedGeneFile.write(json.dumps(updatedGeneDict))

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

def merge_alignments(current_alignment_file, previous_alignment_files, prev_dir):
    # remove extension and split into individual gene names so alignments can be joined file by file
    current_panaroo_label = os.path.basename(current_alignment_file).replace(".aln.fas", "").replace(".fasta", "")
    current_splitNames = current_panaroo_label.split("~~~")
    updated_alignment_files = []
    for fileName in previous_alignment_files:
        previous_panaroo_label = os.path.basename(fileName).replace(".aln.fas", "").replace(".fasta", "")
        previous_splitNames = previous_panaroo_label.split("~~~")
        to_merge = False
        for name in previous_splitNames:
            if name in current_splitNames:
                to_merge = True
                updated_alignment_files.append(current_alignment_file)
        if to_merge:
            current_geneName_set = set(current_splitNames)
            joined_set = current_geneName_set | set(previous_splitNames)
            updatedGeneNames = "~~~".join(list(joined_set))
            updatedFileName = os.path.join(prev_dir, "aligned_gene_sequences", updatedGeneNames + ".aln.fas")
            mafft_command = "mafft --quiet --retree 1 --maxiterate 0 --nofft --add "
            mafft_command += current_alignment_file + " "
            mafft_command += fileName + " > "
            mafft_command += updatedFileName + " "
            mafft_command += " && rm -rf " + fileName # remove outdated alignment file
            subprocess.run(mafft_command, shell=True, check=True)
            return current_alignment_file
        else:
            return ""

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(os.path.join(args.prev_run, os.path.basename(args.accessionFile))):
        sys.stderr.write("\nCopying current run data into " + args.prev_run + "\n")
        subprocess_command = "mkdir " + args.prev_run + " && "
        subprocess_command += "cp -r "
        subprocess_command += args.ncbi_metadata + " "
        subprocess_command += args.geneMetadataDir + " "
        subprocess_command += args.accessionFile + " "
        subprocess_command += args.alignment_dir + " "
        subprocess_command += args.graph_dir + " "
        subprocess_command += args.prev_run
        subprocess.run(subprocess_command, shell=True, check=True)
    else:
        sys.stderr.write("\nOverwriting Panaroo output in " + args.prev_run + " with current merged output\n")
        subprocess.run("cp -rf " + args.graph_dir + " " + args.prev_run, shell=True, check=True)
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
        # copy panarooPairs into previous run dir
        panarooPairs = os.path.join(args.geneMetadataDir, "panarooPairs.json")
        subprocess.run("cp -rf " + panarooPairs + " " + os.path.join(args.prev_run, panarooPairs),shell=True, check=True)
        mergeGeneMetadata(args.geneMetadataDir,
                          args.prev_run)
        sys.stderr.write("\nDone\n")
        # merge mafft alignments of Panaroo output for current and previous runs
        sys.stderr.write("\nMerging current and previous MAFFT alignments\n")
        # parallelise alignment merge
        current_alignment_files = glob.glob(os.path.join(args.alignment_dir, "*.aln.fas")) + glob.glob(os.path.join(args.alignment_dir, "*.fasta"))
        previous_alignment_files = glob.glob(os.path.join(args.prev_run, "aligned_gene_sequences", "*.aln.fas")) + glob.glob(os.path.join(args.prev_run, "aligned_gene_sequences","*.fasta"))
        # import panarooPairs file as there is often an error due to filenames being too long
        with open(os.path.join(geneMetadataDir, "panarooPairs.json"), "r") as jsonFile:
            pairString = jsonFile.read()
        pairs = json.loads(pairString)
        previous_alignment_files_cleaned = []
        for prev_align in previous_alignment_files:
            if not "group_" in prev_align:
                previous_alignment_files_cleaned.append(prev_align)
        current_alignment_subsets = [
            current_alignment_files[i:i + args.n_cpu] for i in range(0, len(current_alignment_files), args.n_cpu)
        ]
        updated_alignment_files = []
        for alignmentFile in tqdm(current_alignment_subsets):
            updated_alignment_files += Parallel(n_jobs=args.n_cpu)(delayed(merge_alignments)(aln,
                                                                                             previous_alignment_files_cleaned,
                                                                                             args.prev_run) for aln in alignmentFile)
        # if current alignment files have not been merged with previous output, copy into prev_dir
        updated_alignment_files = [alnFile for sublist in updated_alignment_files for alnFile in sublist]
        for currentFile in current_alignment_files:
            if not currentFile in updated_alignment_files:
                copy_command = "cp " + currentFile + " " + os.path.join(args.prev_run, "aligned_gene_sequences")
        # genes names with a group prefix change with every run, need to re-align everytime
        sys.stderr.write("\nGenerating MAFFT alignments for unnamed genes\n")
        align_command = "python generate_alignments-runner.py --graph-dir "
        align_command += os.path.join(args.prev_run, "panaroo_output")
        align_command += " --output-dir "
        align_command += os.path.join(args.prev_run, "aligned_gene_sequences")
        align_command += " --post-merge --threads " + str(args.n_cpu)
        subprocess.run(align_command, check=True, shell=True)
        sys.stderr.write("\nDone\n")
        # merge accession IDs for current and previous runs
        sys.stderr.write("\nMerging current and previous accession IDs\n")
        mergeAccessionIDs(args.accessionFile,
                          args.prev_run)
        sys.stderr.write("\nDone\n")
    sys.exit(0)
