import glob
from joblib import Parallel, delayed
import json
import os
import re
from tqdm import tqdm

def get_options():

    import argparse

    description = 'Clean and reformat annotation and sequence files for Panaroo input'
    parser = argparse.ArgumentParser(description=description,
                                        prog='panaroo_clean_inputs')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("-a",
                        "--annotation-dir",
                        dest="annotation_dir",
                        required=True,
                        help='directory of functonal annotations in GFF format for panaroo input',
                        type=str)
    io_opts.add_argument("-g",
                        "--genome-dir",
                        dest="genome_dir",
                        required=True,
                        help='directory of source genomic sequences for annotation files',
                        type=str)
    io_opts.add_argument("-p",
                        "--prodigal-dir",
                        dest="prodigal_dir",
                        required=True,
                        help='directory of prodigal-predicted annotation files',
                        type=str)
    io_opts.add_argument("-i",
                        "--index-file",
                        dest="index_file",
                        required=True,
                        help="JSON file containing integer value to start index from",
                        type=str)
    io_opts.add_argument("-o",
                        "--output-dir",
                        dest="output_dir",
                        required=True,
                        help='directory to output reformatted GFF files',
                        type=str)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads",
                        default=1,
                        type=int)
    args = parser.parse_args()
    return (args)

def translate(seq): #Translate the exon sequence of the gene into its respective amino acid codes using a dictionary
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein =""
    try:
        if len(seq)%3 == 0:
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3]
                protein+= table[codon]
    except:
        protein = ""
    return protein

def reverse_complement(dna):
    try:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        reversed = ''.join([complement[base] for base in dna[::-1]])
    except:
        reversed = ""
    return reversed

def reformat_existing_annotations(stored_annotation,
                                  stored_genome,
                                  COG_no):
    """Reformat existing annotations for compatibility with panaroo"""
    all_region_names = []
    all_region_annotations = []
    for region in range(len(stored_annotation)):
        gene_list = stored_annotation[region].splitlines()
        region_title = gene_list[0].split("\n")[0]
        all_region_names.append("##sequence-region " + region_title)
        for line in range(len(gene_list)):
            if '###' in gene_list[line]:
                gene_list = gene_list[:line]
        gene_list = gene_list[2:]

        genes = []
        for gene_elem in range(len(gene_list)):
            split = re.split(r'\t+', gene_list[gene_elem])
            end = int(split[4])
            start = int(split[3])
            # process_pokka input only looks for CDSs and returns duplicate error when the exon is split
            if split[2] == 'CDS' and (((end - start) + 1) % 3) == 0:
                # add name for consistency if gene is unnamed
                if not ";gene=" in split[8]:
                    gene_list[gene_elem] = gene_list[gene_elem] + ";gene=COG_" + str(COG_no) + ";"
                    COG_no += 1
                genes.append(gene_list[gene_elem])
        # ensure gene ids are unique
        Ids= []
        for gene in range(len(genes)):
            Id = re.search('ID=(.*?);', genes[gene]).group(1)
            Ids.append(Id)
        output_cds = []
        start_index = []
        end_index = []
        sense = []
        for Id in range(len(Ids)):
            if Ids.count(Ids[Id]) == 1:
                output_cds.append(genes[Id])
                no_duplicates_split = re.split(r'\t+', genes[Id])
                start_index.append(int(no_duplicates_split[3]) - 1)
                end_index.append(int(no_duplicates_split[4]) - 1)
                sense.append(no_duplicates_split[6])
        # get genomic sequence corresponding to annotation region
        nucleotides = ''.join(stored_genome[region].split('\n')[1:])
        # remove annotations containing premature stop codons
        to_remove = []
        for sign in range(len(sense)):
            if sense[sign] == '+':
                start_codon = nucleotides[start_index[sign]: start_index[sign] + 3]
                protein = translate(nucleotides[start_index[sign]: end_index[sign] - 2])
                if "*" in protein or protein == "":
                    to_remove.append(sign)
            else:
                reverse_start_codon = nucleotides[end_index[sign] - 2: end_index[sign] + 1]
                protein_reverse = translate(reverse_complement(nucleotides[start_index[sign] + 3: end_index[sign] + 1]))
                if "*" in protein_reverse or protein_reverse == "":
                    to_remove.append(sign)
        for index in sorted(to_remove, reverse=True):
            del output_cds[index]
        if len(output_cds) == 0:
            continue
        cleaned_annotation = "\n".join(str(n) for n in output_cds)
        all_region_annotations.append(cleaned_annotation)
    return all_region_names, all_region_annotations, COG_no

def reformat_predicted_annotations(predicted_annotations, pred_no):
    """Extract region annotation information for prodigal-predicted annotations"""
    predicted_annotations = predicted_annotations.split("# Sequence Data:")
    predicted_region_names = []
    all_predicted_annotations = []
    for region in predicted_annotations:
        if not "##gff-version" in region:
            region_split = region.splitlines()
            predicted_region_annotations = []
            for line in region_split[1:]:
                if not ("# Model Data:" in line or line == "" or line == " "):
                    predicted_region_annotations.append(line + "product=hypothetical protein;gene=PRED_" + str(pred_no) + ";")
                    pred_no += 1
            if not predicted_region_annotations == []:
                all_predicted_annotations.append("\n".join(predicted_region_annotations))
                predicted_region_names.append(line.split("\t")[0])
    return predicted_region_names, all_predicted_annotations, pred_no

def merge_all_annotations(all_region_names,
                          all_region_annotations,
                          predicted_region_names,
                          all_predicted_annotations):
    """Merge prodigal predicted and existing annotations"""
    all_merged_annotations = []
    predicted_regions_merged = []
    # adds predicted annotations to regions that already have existing annotations
    for region_annotation in all_region_annotations:
        region_label = region_annotation.split("\t")[0]
        positions = []
        for existing_line in region_annotation.splitlines():
            tab_split = existing_line.split("\t")
            positions.append((int(tab_split[3]), int(tab_split[4])))
        if region_label in predicted_region_names:
            predicted_region_index = predicted_region_names.index(region_label)
            predicted_annotations = all_predicted_annotations[predicted_region_index]
            annotations_to_add = []
            for predicted_line in predicted_annotations.splitlines():
                predicted_split = predicted_line.split("\t")
                position_tuple = (int(predicted_split[3]), int(predicted_split[4]))
                if not position_tuple in positions:
                    annotations_to_add.append(predicted_line)
            merged_annotation = region_annotation.splitlines() + annotations_to_add
            predicted_regions_merged.append(region_label)
            all_merged_annotations += merged_annotation
    # adds predicted annotations for regions that don't already have existing annotations
    for predicted_region_label in range(len(predicted_region_names)):
        if not predicted_region_names[predicted_region_label] in predicted_regions_merged:
            all_merged_annotations += all_predicted_annotations[predicted_region_label].splitlines()
    return all_region_names, all_merged_annotations

def concatenate_inputs(annotation_file,
                       genome_dir,
                       output_dir,
                       prodigal_dir,
                       index_file):
    """This function merges existing annotations with prodigal-predicted annotations and reformats them for panaroo input"""
    # import annotation and genome files
    label = ('.').join(os.path.basename(annotation_file).split('.')[:-1])
    genome_file = os.path.join(genome_dir, label + ".fna")
    # index used to add names to unnamed genes
    with open(index_file, "r") as indexFile:
        indexNoDict = json.loads(indexFile.read())
    COG_no = int(indexNoDict["cogIndexNo"])
    pred_no = int(indexNoDict["predictedIndexNo"])
    with open(annotation_file, "r") as a:
        stored_annotation = a.read()
    with open(genome_file, "r") as g:
        stored_genome = g.read()
    # if prodigal file exists, import it
    predictedAnnotationFile = os.path.join(prodigal_dir, os.path.basename(annotation_file))
    if os.path.exists(predictedAnnotationFile):
        add_predicted = True
        with open(predictedAnnotationFile, "r") as prodigalFile:
            predicted_annotations = prodigalFile.read()
    else:
        add_predicted = False
    # split annotation and genome into regions
    stored_genome = stored_genome.split('>')[1:]
    stored_annotation = stored_annotation.split("##sequence-region ")[1:]
    all_region_names, all_region_annotations, COG_no = reformat_existing_annotations(stored_annotation,
                                                                                     stored_genome,
                                                                                     COG_no)
    if add_predicted:
        # if predicted annotations exist, reformat them for compatibility
        predicted_region_names, all_predicted_annotations, pred_no = reformat_predicted_annotations(predicted_annotations, pred_no)
        # if predicted annotations exist, merge them with the existing annotations
        all_region_names, all_region_annotations = merge_all_annotations(all_region_names,
                                                                         all_region_annotations,
                                                                         predicted_region_names,
                                                                         all_predicted_annotations)
    # write out reformatted GFF file
    annotated_file = "\n".join(all_region_names + all_region_annotations) + "\n##FASTA\n" + ">".join(stored_genome)
    filename_cleaned = os.path.join(output_dir, os.path.basename(annotation_file))
    with open(filename_cleaned,'w') as o:
        o.write(annotated_file)
    # update cog index number for subsequent runs
    indexNoDict["cogIndexNo"] = COG_no
    indexNoDict["predictedIndexNo"] = pred_no
    with open(index_file, "w") as indexFile:
        indexFile.write(json.dumps(indexNoDict))

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()
    # create output directory
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    annotation_files = glob.glob(os.path.join(args.annotation_dir, "*.gff"))
    # parallelise
    job_list = [
        annotation_files[i:i + args.n_cpu] for i in range(0, len(annotation_files), args.n_cpu)
    ]
    # iterate through annotations and reformat for panaroo input
    for job in tqdm(job_list):
        Parallel(n_jobs=args.n_cpu)(delayed(concatenate_inputs)(annotation,
                                                                args.genome_dir,
                                                                args.output_dir,
                                                                args.prodigal_dir,
                                                                args.index_file) for annotation in job)
if __name__ == '__main__':
    main()
