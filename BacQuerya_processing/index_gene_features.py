"""
Build a searcheable COBS index from the output of feature_extract
"""
import cobs_index as cobs
import glob
from joblib import Parallel, delayed
import json
import os
import shutil
import sys
from tqdm import tqdm
import tempfile

def get_options():

    import argparse
    description = 'Generate a searcheable COBS index from genetic features'
    parser = argparse.ArgumentParser(description=description,
                                     prog='index-genes')
    io_opts = parser.add_argument_group('Inputs')

    io_opts.add_argument("-i",
                        "--input-file",
                        dest="input_file",
                        required=True,
                        help='file of all genetic sequences output by feature_extract',
                        type=str)
    io_opts.add_argument("-o",
                        "--output-dir",
                        dest="output_dir",
                        required=True,
                        help="output directory for generated index",
                        type=str)
    io_opts.add_argument("--kmer-length",
                        dest="kmer_length",
                        help="specify kmer length for the constructed index (default = 15)",
                        default=15,
                        type=int)
    io_opts.add_argument("--false-positive",
                        dest="fpr",
                        help="false positive rate for index. Greater fpr means smaller index (default = 0.01).",
                        default=0.01,
                        type=float)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    io_opts.add_argument("--dirty",
                        dest="dirty",
                        help="keep gene files used to build the index",
                        action='store_true',
                        default=False)

    args = parser.parse_args()
    return (args)

def write_files(gene_dict, temp_dir):
    """Write gene sequences to individual files with index as filename"""
    gene_index = str(gene_dict["index"])
    gene_sequence = gene_dict["sequence"]
    with open(os.path.join(temp_dir, gene_index + ".txt"), "w") as g:
        g.write(gene_sequence)

def create_index(temp_dir, output_dir, kmer_length, fpr):
    """"Create COBS index"""
    params = cobs.CompactIndexParameters()
    params.term_size = kmer_length
    params.clobber = True               # overwrite output and temporary files
    params.false_positive_rate = fpr    # higher false positive rate -> smaller index
    cobs.compact_construct(os.path.join(temp_dir),
                           os.path.join(output_dir,
                                        str(kmer_length) + "_index.cobs_compact"),
                           index_params=params)

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
    with open(args.input_file, "r") as f:
        gene_dicts_str = f.read()
    gene_dicts = json.loads(gene_dicts_str)
    gene_dicts = gene_dicts["information"]
    job_list = [
        gene_dicts[i:i + args.n_cpu] for i in range(0, len(gene_dicts), args.n_cpu)
    ]
    # parrallelise writing of gene-specific files for indexing
    for job in tqdm(job_list):
        Parallel(n_jobs=args.n_cpu)(delayed(write_files)(g,
                                                        temp_dir) for g in job)
    create_index(temp_dir,
                 args.output_dir,
                 args.kmer_length,
                 args.fpr)
    if not args.dirty:
        shutil.rmtree(temp_dir)
    sys.exit(0)

if __name__ == '__main__':
    main()
