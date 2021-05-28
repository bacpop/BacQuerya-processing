import json
import sys
from tqdm import tqdm

def get_options():

    import argparse

    description = 'Download reads of interest by accession ID'
    parser = argparse.ArgumentParser(description=description,
                                        prog='extract_entrez_reads')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("--file",
                        dest="input",
                        required=True,
                        help='txt file listing isolate biosamples and EBI assembly paths (http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/)',
                        type=str)
    args = parser.parse_args()
    return (args)

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()
    with open(args.input, "r") as inFile:
        assembly_list = inFile.read().splitlines()
    assembly_link_dict = {}
    for line in tqdm(assembly_list):
        split = line.split("\t")
        biosample = split[0]
        assemblyURL = split[1].replace("/ebi/ftp/", "http://ftp.ebi.ac.uk/")
        assembly_link_dict[biosample] = assemblyURL
    with open("661K_biosampleAssemblyURLs.json", "w") as keyFile:
        keyFile.write(json.dumps(assembly_link_dict))
    sys.exit(0)

if __name__ == '__main__':
    main()
