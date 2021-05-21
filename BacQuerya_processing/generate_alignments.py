import glob
from joblib import Parallel, delayed
import json
import networkx as nx
import os
import pandas as pd
import shutil
import subprocess
import sys
from tqdm import tqdm

def get_options():
    import argparse
    description = 'Extract features from gff and sequence files'
    parser = argparse.ArgumentParser(description=description,
                                        prog='extract_genes')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-g",
                        "--graph-dir",
                        dest="graph_dir",
                        required=True,
                        help="directory of Panaroo graph",
                        type=str)
    io_opts.add_argument("-e",
                        "--extracted-genes",
                        dest="extracted_genes",
                        required=True,
                        help="directory output by extract_genes.py for current run",
                        type=str)
    io_opts.add_argument("-o",
                        "--output-dir",
                        dest="output_dir",
                        required=True,
                        help="output directory",
                        type=str)
    io_opts.add_argument("--post-merge",
                        dest="post_merge",
                        help="alignement script is being run after panaroo_merge",
                        action='store_true',
                        default=False)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    args = parser.parse_args()
    return (args)

def mafft_align(input_file, output_dir):
    """Multi-threaded mafft alignment using the subprocess module"""
    subprocess_command = "mafft --quiet --retree 1 --maxiterate 0 --nofft --parttree "
    subprocess_command += input_file + " > "
    subprocess_command += os.path.join(output_dir, os.path.basename(input_file).replace(".fasta", ".aln.fas"))
    subprocess.run(subprocess_command, shell=True, check=True)

def deleteFile(filename):
    if os.path.isdir(filename):
        shutil.rmtree(filename)

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    raw_files = os.path.join(args.output_dir, "raw_files")
    if not os.path.exists(raw_files):
        os.mkdir(raw_files)
    gene_data = pd.read_csv(os.path.join(args.graph_dir, "gene_data.csv"))
    # import panarooPairs file as there is often an error due to filenames being too long
    with open(os.path.join(args.extracted_genes, "panarooPairs.json"), "r") as jsonFile:
        pairString = jsonFile.read()
    pairs = json.loads(pairString)
    # convert gene data file to dict to improve efficiency
    sys.stderr.write("\nConverting gene_data.csv to dict\n")
    gene_data_json = {}
    for row in tqdm(range(len(gene_data["clustering_id"]))):
        sequence = gene_data["dna_sequence"][row]
        cluster_dict = {gene_data["clustering_id"][row] : sequence}
        gene_data_json.update(cluster_dict)
    # write msa file for mafft input
    sys.stderr.write('\nWriting unaligned multi-fasta files from Panaroo output\n')
    G = nx.read_gml(os.path.join(args.graph_dir, "final_graph.gml"))
    unaligned_files = []
    if not args.post_merge:
        for node in tqdm(G._node):
            y = G._node[node]
            gene_names = y["name"]
            # unnamed genes are assigned non-consistent names with a group_ prefix
            if not "group_" in y["name"]:
                for key, value in pairs.items():
                    if gene_names == value["panarooNames"]:
                        consistent_name = value["consistentNames"]
                if not len(y["members"]) == 1:
                    multiFSAline = []
                    for mem in range(len(y["members"])):
                        isol_label = G.graph["isolateNames"][mem]
                        member_sequence = gene_data_json[y["geneIDs"].split(";")[mem]]
                        multiFSAline.append(">" + isol_label + "\n" + member_sequence)
                    # use consistentName as filename
                    filename = os.path.join(raw_files, consistent_name + ".fasta")
                    unaligned_files.append(filename)
                    with open(filename, "w") as outFile:
                        outFile.write("\n".join(multiFSAline))
                else:
                    isol_label = G.graph["isolateNames"][y["members"][0]]
                    member_sequence = gene_data_json[y["geneIDs"].split(";")[0]] # just in case we are looking at a paralog
                    multiFSAline = ">" + isol_label + "\n" + member_sequence
                    with open(os.path.join(args.output_dir, gene_names + ".fasta"), "w") as outFile:
                        outFile.write(multiFSAline)
    else:
        # delete group alignments if they are already present- the gene names are incorrect
        group_files = glob.glob(os.path.join(args.output_dir, "group_*"))
        for groupFile in group_files:
            deleteFile(groupFile)
        # iterate through group nodes
        for node in tqdm(G._node):
            y = G._node[node]
            gene_names = y["name"]
            # unnamed genes are assigned inconsistent names with a group_ prefix.
            if "group_" in y["name"]:
                if not len(y["members"]) == 1:
                    multiFSAline = []
                    for mem in range(len(y["members"])):
                        isol_label = G.graph["isolateNames"][mem]
                        member_sequence = gene_data_json[y["geneIDs"].split(";")[mem]]
                        multiFSAline.append(">" + isol_label + "\n" + member_sequence)
                    filename = os.path.join(raw_files, gene_names + ".fasta")
                    unaligned_files.append(filename)
                    with open(filename, "w") as outFile:
                        outFile.write("\n".join(multiFSAline))
                else:
                    isol_label = G.graph["isolateNames"][y["members"][0]]
                    member_sequence = gene_data_json[y["geneIDs"].split(";")[0]] # just in case we are looking at a paralog
                    multiFSAline = ">" + isol_label + "\n" + member_sequence
                    with open(os.path.join(args.output_dir, gene_names + ".fasta"), "w") as outFile:
                        outFile.write(multiFSAline)
    # run mafft on the unaligned MSA files
    sys.stderr.write('\nAligning gene multi-fasta files\n')
    # parallelise mafft alignment
    file_list = [
        unaligned_files[i:i + args.n_cpu] for i in range(0, len(unaligned_files), args.n_cpu)
    ]
    for fasta_file in tqdm(file_list):
        Parallel(n_jobs=args.n_cpu)(delayed(mafft_align)(in_file,
                                                         args.output_dir) for in_file in fasta_file)
    # clean up unaligned files
    shutil.rmtree(raw_files)
    sys.exit(0)

if __name__ == '__main__':
    main()