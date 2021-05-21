"""Replace the name of downloaded GPS assemblies and assembly metadata with the lane ID"""
import glob
import os
import sys

def get_options():

    import argparse

    description = 'Download reads of interest by accession ID'
    parser = argparse.ArgumentParser(description=description,
                                        prog='extract_entrez_reads')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("--input-dir",
                        dest="input_dir",
                        required=True,
                        help='directory containing files to rename',
                        type=str)
    args = parser.parse_args()
    return (args)

def main():
    args = get_options()
    retieved_files = glob.glob(os.path.join(args.input_dir, "*"))
    for file in retieved_files:
        dirname = os.path.dirname(file)
        print(dirname)
        newFilename = ".".join(file.split(".")[1:]).split("_")[1:]
        newFilename = os.path.join(dirname, "_".join(newFilename[:2]) + "#" + newFilename[2])
        os.rename(file, newFilename)
    sys.exit(0)

if __name__ == '__main__':
    main()