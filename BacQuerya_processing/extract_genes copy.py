#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses GFF and Sequence parsing to convert GFF files to json strings and extract gene names and sequences for COBS indexing.
"""
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from elasticsearch import Elasticsearch
from joblib import Parallel, delayed
import json
import glob
import networkx as nx
import os
import pandas as pd
import re
import requests
import sys
from tqdm import tqdm
import xmltodict

from BacQuerya_processing.index_gene_features import elasticsearch_isolates
from BacQuerya_processing.panaroo_clean_inputs import reverse_complement, translate

def get_options():

    import argparse

    description = 'Extract features from gff and sequence files'
    parser = argparse.ArgumentParser(description=description,
                                        prog='extract_genes')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-s",
                        "--sequences",
                        dest="sequences",
                        required=True,
                        help='directory of genomic sequences',
                        type=str)
    io_opts.add_argument("-a",
                        "--gffs",
                        dest="gffs",
                        required=True,
                        help="directory of annotations in gff format",
                        type=str)
    io_opts.add_argument("-g",
                        "--graph-dir",
                        dest="graph_dir",
                        required=True,
                        help="directory of Panaroo graph",
                        type=str)
    io_opts.add_argument("-k",
                        "--isolate-pairs",
                        dest="isolate_pairs",
                        required=True,
                        help="json file of index isolate name pairs",
                        type=str)
    io_opts.add_argument("-j",
                        "--isolate-json",
                        dest="isolate_json",
                        required=True,
                        help="json of all isolate attributes output by extract_assembly_stats",
                        type=str)
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,
                        help="output directory",
                        type=str)
    io_opts.add_argument("-i",
                        "--index-file",
                        dest="index_file",
                        required=True,
                        help="JSON file containing integer value to start index from",
                        type=str)
    io_opts.add_argument("--elastic-index",
                        dest="elastic",
                        help="don't write gene json and index directly in script",
                        action='store_true',
                        default=False)
    io_opts.add_argument("--index-name",
                        dest="index_name",
                        required=False,
                        help="index to create/append to",
                        type=str)
    io_opts.add_argument("--biosampleJSON",
                        dest="biosampleJSON",
                        required=True,
                        help="json file of biosample isolate name pairs",
                        type=str)
    io_opts.add_argument("--run-type",
                        dest="run_type",
                        required=True,
                        help="whether pipeline is indexing reference or query isolates",
                        type=str)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    args = parser.parse_args()
    return (args)


def searchPfam(proteinSequence):
    """Search for hypothetical proteins in the Pfam database"""
    headers= {"Expect": "", "Accept": "text/xml"}
    parameters = { 'hmmdb' : 'pfam', 'seq': proteinSequence}
    url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'
    urlResponse = requests.post(url, headers = headers, data = parameters)
    try:
        pfamResult = xmltodict.parse(urlResponse.text)
    except:
        urlText = urlResponse.text
        # xml for hits include tags of integers that xmltodict does not recoginise
        urlText = urlText.split("<")
        cleanedText = []
        for tag in urlText:
            if not tag.split(" ")[0].isdigit():
                cleanedText.append(tag)
        xmlCleaned = "<".join(cleanedText)
        pfamResult = xmltodict.parse(xmlCleaned)
    if "hits" in pfamResult["opt"]["data"].keys():
        match = pfamResult["opt"]["data"]["hits"]
        if not isinstance(match, list):
            pfamDict = {"pfam_names": match["@name"],
                        "pfam_accessions": match["@acc"],
                        "pfam_biases": match["@bias"],
                        "pfam_descriptions": match["@desc"],
                        "pfam_evalues": match["@evalue"]}
        else:
            names = []
            accessions = []
            biases = []
            descriptions = []
            evalues = []
            for hit in match:
                names.append(hit["@name"])
                accessions.append(hit["@acc"])
                biases.append(hit["@bias"])
                descriptions.append(hit["@desc"])
                evalues.append(hit["@evalue"])
            pfamDict = {"pfam_names": names,
                        "pfam_accessions": accessions,
                        "pfam_biases": biases,
                        "pfam_descriptions": descriptions,
                        "pfam_evalues": evalues}
        return pfamDict
    return None

def update_panaroo_outputs(G,
                           panarooNames,
                           updatedPanarooNames,
                           updatedDescriptions,
                           graph_dir,
                           gene_data):
    """Use consistentNames/ pfam names and pfam descriptions to update all relevant panaroo outputs.
       This makes it significantly easier to keep track of everything that's going on"""
    sys.stderr.write('\nUpdating Panaroo graph\n')
    # let's start by updating the graph itself
    cluster_name_dict = {}
    for node in tqdm(G._node):
        y = G._node[node]
        gene_names = y["name"]
        splitNames = gene_names.split("~~~")
        if not all("PRED_" in name for name in splitNames):
            node_index = panarooNames.index(gene_names)
            node_gene_name = updatedPanarooNames[node_index]
            node_annotation = updatedDescriptions[node_index]
            # used to update gene_data.csv
            cluster_ids = y["geneIDs"].split(";")
            for clus in cluster_ids:
                cluster_name_dict.update({clus: {"name": node_gene_name, "description": node_annotation}})
            # update node name
            y["name"] = node_gene_name
            y["description"] = node_annotation
    # overwrite the previous graph
    nx.write_gml(G, os.path.join(graph_dir, "final_graph.gml"))
    del G
    # update gene_data.csv
    sys.stderr.write('\nUpdating all gene_data.csv file\n')
    for row in tqdm(range(len(gene_data["clustering_id"]))):
        if gene_data["clustering_id"][row] in cluster_name_dict.keys():
            gene_data["gene_name"][row] = cluster_name_dict[gene_data["clustering_id"][row]]["name"]
            gene_data["description"][row] = cluster_name_dict[gene_data["clustering_id"][row]]["description"]
    gene_data.to_csv(os.path.join(graph_dir, "gene_data.csv"), index=False)
    del gene_data
    # update gene_presence_absence.csv and gene_presence_absence.rtab
    sys.stderr.write('\nUpdating gene_presence_absence files\n')
    gene_presence_absence = pd.read_csv(os.path.join(graph_dir, "gene_presence_absence.csv"))
    roary_presence_absence = pd.read_csv(os.path.join(graph_dir, "gene_presence_absence_roary.csv"))
    updated_cog_dict = {}
    for name in tqdm(range(len(panarooNames))):
        # update csv file
        annotation_row_index = list(gene_presence_absence["Gene"]).index(panarooNames[name])
        non_unique = ";".join(updatedPanarooNames[name].split("~~~"))
        gene_presence_absence["Gene"][annotation_row_index] = updatedPanarooNames[name]
        gene_presence_absence["Non-unique Gene name"][annotation_row_index] = non_unique
        gene_presence_absence["Annotation"][annotation_row_index] = updatedDescriptions[name]
        # update roary csv file
        roary_row_index = list(roary_presence_absence["Gene"]).index(panarooNames[name])
        roary_presence_absence["Gene"][roary_row_index] = updatedPanarooNames[name]
        roary_presence_absence["Non-unique Gene name"][roary_row_index] = non_unique
        roary_presence_absence["Annotation"][roary_row_index] = updatedDescriptions[name]
        updated_cog_dict.update({panarooNames[name]: updatedPanarooNames[name]})
    sys.stderr.write("\nWriting csv files\n")
    gene_presence_absence.to_csv(os.path.join(graph_dir, "gene_presence_absence.csv"), index=False)
    roary_presence_absence.to_csv(os.path.join(graph_dir, "gene_presence_absence_roary.csv"), index=False)
    del gene_presence_absence
    del roary_presence_absence
    sys.stderr.write('\nUpdating gene_presence_absence.Rtab\n')
    with open(os.path.join(graph_dir, "gene_presence_absence.Rtab"), "r") as inRtab:
        rtab_pres_abs = inRtab.read().splitlines()
    for key, value in tqdm(updated_cog_dict.items()):
        # update RTAB file
        for line in range(len(rtab_pres_abs)):
            cog_names = rtab_pres_abs[line].split("\t")[0]
            split_cog_names = cog_names.split("~~~")
            key_split = key.split("~~~")
            if all(k_name in split_cog_names for k_name in key_split):
                rtab_pres_abs[line] = rtab_pres_abs[line].replace(cog_names, value)
    with open(os.path.join(graph_dir, "gene_presence_absence.Rtab"), "w") as outRtab:
        outRtab.write("\n".join(rtab_pres_abs))
    del rtab_pres_abs
    sys.stderr.write('\nUpdating pan_genome_reference.fa\n')
    with open(os.path.join(graph_dir, "pan_genome_reference.fa"), "r") as refFile:
        pan_genome_reference = refFile.read().split(">")
    for key, value in tqdm(updated_cog_dict.items()):
        # update pan_genome_reference.fa
        for line in range(len(pan_genome_reference)):
            cog_names = pan_genome_reference[line].split("\n")[0]
            split_cog_names = cog_names.split("~~~")
            k_split = key.split("~~~")
            if all(k in split_cog_names for k in k_split):
                pan_genome_reference[line] = pan_genome_reference[line].replace(cog_names, value)
    pan_genome_reference = ">".join(pan_genome_reference)
    with open(os.path.join(graph_dir, "pan_genome_reference.fa"), "w") as outRefFile:
        outRefFile.write(pan_genome_reference)
    del pan_genome_reference
    sys.stderr.write('\nUpdating struct_presence_absence.Rtab\n')
    with open(os.path.join(graph_dir, "struct_presence_absence.Rtab"), "r") as inStruct:
        struct_pres_abs = inStruct.read().splitlines()
    # update struct_presence_absence.rtab
    for line in tqdm(range(len(struct_pres_abs))):
        cog_names = struct_pres_abs[line].split("\t")[0]
        split_cog_names = cog_names.split("-")
        for key, value in updated_cog_dict.items():
            for name in split_cog_names:
                if name == key:
                    struct_pres_abs[line] = struct_pres_abs[line].replace(name, value)
    struct_pres_abs = "\n".join(struct_pres_abs)
    with open(os.path.join(graph_dir, "struct_presence_absence.Rtab"), "w") as outStructFile:
        outStructFile.write(struct_pres_abs)
    del struct_pres_abs

def generate_library(graph_dir,
                     index_no,
                     output_dir,
                     isolateIndexJSON,
                     threads,
                     biosampleJSON):
    """Extract all newly annotated/identified genes from panaroo graph and update geneJSON"""
    G = nx.read_gml(os.path.join(graph_dir, "final_graph.gml"))
    gene_data = pd.read_csv(os.path.join(graph_dir, "gene_data.csv"))
    # convert gene_data df to json to speed up sequence extraction
    sys.stderr.write('\nConverting gene data dataframe to JSON\n')
    gene_data_json = {}
    updated_genes = []
    for row in tqdm(range(len(gene_data["clustering_id"]))):
        sequence = gene_data["dna_sequence"][row]
        cluster_dict = {gene_data["clustering_id"][row] : gene_data["annotation_id"][row]}
        gene_data_json.update(cluster_dict)
    with open(biosampleJSON, "r") as bios:
        label_accession_str = bios.read()
    label_accession_pairs = json.loads(label_accession_str)
    num_isolates = len(G.graph["isolateNames"])
    annotationID_key_updated_genes = {}
    # need to update panarooPairs.json if it already exists, if not then create it.
    previousRunPanarooPairs = os.path.join(os.path.dirname(graph_dir), os.path.basename(output_dir), "panarooPairs.json")
    if os.path.exists(previousRunPanarooPairs):
        with open(previousRunPanarooPairs, "r") as prevFile:
            panaroo_pairsJSON = prevFile.read()
        panaroo_pairs = json.loads(panaroo_pairsJSON)
        update_index_no = False
    else:
        panaroo_pairs = {}
        update_index_no = True
    # iterate through panaroo graph to extract gene information if node is not present in panarooPairs or has been updated
    sys.stderr.write('\nExtracting node information from Panaroo graph\n')
    updatedPanarooNames = []
    currentPanarooNames = []
    updatedDescriptions = []
    all_names_set = set()
    for node in tqdm(G._node):
        y = G._node[node]
        gene_names = y["name"]
        splitNames = gene_names.split("~~~")
        frequency = round((len(y["members"])/num_isolates)*100, 1)
        if not update_index_no:
            # this checks if any of the nodes gene names are in the panaroo pair values
            # if so the index no will be the key
            # if not the index no will be the current index no in index_values
            for key, value in panaroo_pairs.items():
                if any(name in value["panarooNames"] for name in splitNames):
                    index_no = int(key)
                    all_names_set.add(value["consistentNames"])
                else:
                    index_no = index_no
        if not all("PRED_" in name for name in splitNames):
            member_labels = []
            annotation_ids = []
            biosample_labels = []
            for mem in range(len(y["members"])):
                isol_label = G.graph["isolateNames"][mem]
                member_labels.append(isol_label)
                biosample_labels.append(label_accession_pairs[isol_label])
                member_annotation_id = gene_data_json[y["geneIDs"].split(";")[mem]]
                annotation_ids.append(member_annotation_id)
            isolate_indices = [isolateIndexJSON[label] for label in member_labels]
            # supplement annotation with pfam search result. Tend to be more up to date
            if not (y["description"] == "" or y["description"] == "hypothetical protein" or y["description"] == "Hypothetical protein"):
                panarooDescription = y["description"].split(";")
                updatedDescriptions.append(";".join(panarooDescription).replace(",", ";"))
                pfamResult = None
            else:
                panarooDescription = ["Hypothetical protein"]
                pfamResult = searchPfam(y["protein"].split(";")[0])
                if pfamResult:
                    # these are new descriptions identified by searching pfam
                    pfamDescriptions = pfamResult["pfam_descriptions"]
                    if isinstance(pfamDescriptions, str):
                        pfamDescriptions = [pfamDescriptions]
                    updatedDescriptions.append(";".join(panarooDescription + pfamDescriptions).replace(",", ";"))
                else:
                    updatedDescriptions += panarooDescription
            # this is the dictionary used to inform the geneDisplay page for the frontend
            annotation_dict = {"panarooDescriptions" : panarooDescription,
                               "panarooFrequency": frequency,
                               "gene_index": index_no,
                               "foundIn_labels": member_labels,
                               "foundIn_indices": isolate_indices,
                               "foundIn_biosamples": biosample_labels,
                               "member_annotation_ids": annotation_ids}
            if not update_index_no:
                consistent_name = panaroo_pairs[index_no]["consistentNames"]
            else:
                consistent_name = "COG_" + str(index_no)
            if pfamResult:
                annotation_dict.update(pfamResult)
            reject_list = ["group_", "PRED_"]
            if not (all(any(x in name for x in reject_list) for name in splitNames)):
                # sets the consistent name to one of the splitnames
                newSplitNames = []
                for name in splitNames:
                    if not ("PRED_" in name or name in all_names_set) and update_index_no:
                        consistent_name = name
                        newSplitNames.append(name)
                        all_names_set.add(name)
                if newSplitNames == []:
                    newSplitNames = [consistent_name]
                # we have removed PRED_, UNNAMED_ and group_ from the name for subsequent panaroo runs
                newGeneNames = "~~~".join(newSplitNames)
            else:
                # apply a consistent name if all annotations are named with an UNNAMED_ prefix
                newGeneNames = consistent_name
            panaroo_pairs.update({index_no: {"panarooNames": newGeneNames, "consistentNames": consistent_name}})
            currentPanarooNames.append(gene_names)
            # we are going to update the node names in the graph
            updatedPanarooNames.append(newGeneNames)
            # prevent duplicate names being returned from gene search results
            all_names_set.add(consistent_name)
            annotation_dict.update({"panarooNames": newGeneNames,"consistentNames": consistent_name})
            # add annotation dict to a list that is then saved. This is as a backup if our indexed data is lost and can be directly indexed by index_gene_features.py.
            updated_genes.append(annotation_dict)
            for annot_ID in annotation_ids:
                annotationID_key_updated_genes.update({annot_ID: annotation_dict})
            index_no += 1
    # update all of the panaroo outputs with information sourced above
    update_panaroo_outputs(G,
                           currentPanarooNames,
                           updatedPanarooNames,
                           updatedDescriptions,
                           graph_dir,
                           gene_data)
    # write name, index pairs in graph for COBS indexing in index_gene_features
    with open(os.path.join(output_dir, "panarooPairs.json"), "w") as o:
        o.write(json.dumps(panaroo_pairs))
    return annotationID_key_updated_genes, updated_genes, index_no

def build_gff_jsons(gff_file,
                    seq_dir,
                    isolateIndexJSON):
    """Convert isolate GFF files to JSON strings and identify the corresponding genomic sequence"""
    label = os.path.basename(gff_file.replace(".gff", ""))
    in_seq_file = os.path.join(seq_dir, label + ".fna")
    in_seq_handle = open(in_seq_file,'r')
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()
    json_feature_dicts = {}
    json_feature_list = []
    with open(gff_file, "r") as g:
        gff_content = g.read().split("##")[2:-1]
    gff_region_title = ""
    for region in gff_content:
        if "sequence-region" == region.split(" ")[0]:
            gff_region_title = region.split(" ")[1]
        region_annotations = region.splitlines()
        for annotation_line in region_annotations:
            if gff_region_title == annotation_line.split("\t")[0]:
                region_sequence = seq_dict[gff_region_title]
                annotation_content = annotation_line.split("\t")
                qualifiers = annotation_content[8]
                start = int(annotation_content[3])
                end = int(annotation_content[4])
                strand = annotation_content[6]
                if strand == "+":
                    feature_sequence = str(region_sequence[start : end])
                    strand = 1
                elif strand == "-":
                    feature_sequence = str(region_sequence[start : end].reverse_complement())
                    strand = -1
                isolate_index = isolateIndexJSON[label]
                json_features = {"type":annotation_content[2],
                                "strand":strand,
                                "start":start,
                                "end":end,
                                "sequenceLength":len(feature_sequence),
                                "isolateName": label,
                                "isolateIndex": isolate_index}
                qualifiers = qualifiers.split(";")
                for qual in qualifiers:
                    split = qual.split("=")
                    attribute = split[0]
                    val = split[1]
                    json_features.update({attribute : val})
                json_feature_list.append(json_features)
    json_feature_dicts.update({label : json_feature_list})
    return json_feature_dicts

def append_gene_indices(isolate_file, all_features):
    """Add indices of all genes within the isolate to the isolate attribute json"""
    with open(isolate_file, "r") as f:
        isolate_json = f.read()
    isolate_dict = json.loads(isolate_json)
    for isol_name in tqdm(range(len(isolate_dict["information"]))):
        isolate_gene_indices = []
        isolate_gene_names = []
        isolate_non_CDS = []
        isolateMetadataName = isolate_dict["information"][isol_name]["isolateNameUnderscore"]
        for annotation_line in all_features:
            if annotation_line["isolateName"] == isolateMetadataName:
                if "panarooNames" in annotation_line.keys():
                    isolate_gene_indices.append(annotation_line["gene_index"])
                    isolate_gene_names.append(annotation_line["consistentNames"])
                if "featureIndex" in annotation_line.keys():
                    isolate_non_CDS.append(annotation_line["featureIndex"])
        if not len(isolate_gene_names) == 0:
            #isolate_dict["information"][isol_name]["geneIndices"] = isolate_gene_indices
            isolate_dict["information"][isol_name]["consistentNames"] = isolate_gene_names
            #isolate_dict["information"][isol_name]["nonCDSIndices"] = isolate_non_CDS
    with open(isolate_file, "w") as o:
        o.write(json.dumps(isolate_dict))

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # load isolate-index k, v pairs
    with open(args.isolate_pairs, "r") as isolateIndex:
        isolateString = isolateIndex.read()
    isolateIndexJSON = json.loads(isolateString)
    # standardise annotations from panaroo output
    with open(args.index_file, "r") as indexFile:
        indexNoDict = json.loads(indexFile.read())
    index_no = int(indexNoDict["geneIndexNo"])
    sys.stderr.write('\nLoading Panaroo graph\n')
    annotationID_key_updated_genes, updated_annotations, index_no = generate_library(args.graph_dir,
                                                                                     index_no,
                                                                                     args.output_dir,
                                                                                     isolateIndexJSON,
                                                                                     args.n_cpu,
                                                                                     args.biosampleJSON)
    gffs = glob.glob(args.gffs + '/*.gff')
    sys.stderr.write('\nConverting annotation files to JSON\n')
    gff_list = [
        gffs[i:i + args.n_cpu] for i in range(0, len(gffs), args.n_cpu)
    ]
    # parrallelise conversion of GFF to json
    all_genes = []
    json_feature_dicts = {}
    for gff in tqdm(gff_list):
        features = Parallel(n_jobs=args.n_cpu)(delayed(build_gff_jsons)(g,
                                                                        args.sequences,
                                                                        isolateIndexJSON) for g in gff)
        for feature_item in features:
            json_feature_dicts.update(feature_item)
    sys.stderr.write('\nUpdating annotation JSON with Panaroo-sourced information\n')
    all_features = []
    for isolate_label, annotations in tqdm(json_feature_dicts.items()):
        non_CDS_index_no = index_no
        for annotation_line in annotations:
            if annotation_line["type"] == "CDS" and annotation_line["ID"] in annotationID_key_updated_genes.keys():
                annotation_line.update(annotationID_key_updated_genes[annotation_line["ID"]])
                all_features.append(annotation_line)
            else:
                annotation_line.update({"featureIndex": str(annotation_line["isolateIndex"]) + "_" + str(non_CDS_index_no)})
                all_features.append(annotation_line)
                non_CDS_index_no += 1
    # add gene indices to isolate jsons
    sys.stderr.write('\nAdding gene indices to isolate assembly JSONs\n')
    append_gene_indices(args.isolate_json,
                        all_features)
    if args.elastic:
        # directly add information to elasticindex
        sys.stderr.write('\nBuilding Elastic Search index\n')
        #elasticsearch_isolates(updated_annotations, args.index_name)
    sys.stderr.write('\nWriting gene JSON files\n')
    with open(os.path.join(args.output_dir, "annotatedNodes.json"), "w") as n:
        n.write(json.dumps({"information":updated_annotations}))
    # update isolate index number for subsequent runs
    indexNoDict["geneIndexNo"] = index_no
    with open(args.index_file, "w") as indexFile:
        indexFile.write(json.dumps(indexNoDict))
    sys.exit(0)

if __name__ == '__main__':
    main()
