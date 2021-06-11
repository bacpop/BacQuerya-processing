import json
import sys
from tqdm import tqdm

def get_options():

    import argparse

    parser = argparse.ArgumentParser(prog='extract_entrez_reads')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("--metadata-file",
                        dest="metadata",
                        required=True,
                        help='Isolate asssembly metadata output by extact_asembly_stats.py',
                        type=str)
    args = parser.parse_args()
    return (args)

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()
    sys.stderr.write("\nImporting metadata file\n")
    with open(args.metadata, "r") as metaFile:
        metadata = json.loads(metaFile.read())["information"]
    sys.stderr.write("\nIterating through isolates and extracting assembly stats\n")
    data = {}
    for isolate in tqdm(metadata):
        isolate_stats = {}
        isolate_stats["contig_N50"] = isolate["contig_stats"]["N50"]
        isolate_stats["gc_content"] = isolate["contig_stats"]["gc_content"]
        isolate_stats["contig_mean_length"] = isolate["contig_stats"]["mean"]
        isolate_stats["contig_count"] = isolate["contig_stats"]["sequence_count"]
        isolate_stats["genome_length"] = isolate["contig_stats"]["total_bps"]
        isolate_stats["scaffold_N50"] = isolate["contig_stats"]["N50"]
        isolate_stats["scaffold_mean_length"] = isolate["scaffold_stats"]["mean"]
        isolate_stats["scaffold_count"] = isolate["scaffold_stats"]["sequence_count"]
        data[isolate["isolateNameUnderscore"]] = isolate_stats
    with open("population_assembly_stats.json", "w") as outJSON:
        outJSON.write(json.dumps(data))
    sys.exit(0)

if __name__ == '__main__':
    main()