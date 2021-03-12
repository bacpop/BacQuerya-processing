import glob
from joblib import Parallel, delayed
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

def concatenate_inputs(annotation_file, genome_dir, output_dir):
    """Most annotations consist of multiple contigs and Panaroo only accepts prokka-formatted genomic regions"""
    # import annotation and genome files
    label = ('.').join(os.path.basename(annotation_file).split('.')[:-1])
    genome_file = os.path.join(genome_dir, label + ".fna")
    with open(annotation_file, "r") as a:
        stored_annotation = a.read()
    with open(genome_file, "r") as g:
        stored_genome = g.read()
    # split annotation and genome into regions
    stored_genome = stored_genome.split('>')[1:]
    stored_annotation = stored_annotation.split("##sequence-region ")[1:]
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
    # write out reformatted GFF file
    annotated_file = "\n".join(all_region_names + all_region_annotations) + "\n##FASTA\n" + ">".join(stored_genome)
    filename_cleaned = os.path.join(output_dir, os.path.basename(annotation_file))
    with open(filename_cleaned,'w') as o:
        o.write(annotated_file)
    return

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
                                                                args.output_dir) for annotation in job)
if __name__ == '__main__':
    main()
