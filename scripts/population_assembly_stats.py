import json
import sys
from tqdm import tqdm

def get_options():

    import argparse

    parser = argparse.ArgumentParser(prog='extract_entrez_reads')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("--assembly-metadata",
                        dest="assemblies",
                        required=True,
                        help='Isolate assembly metadata output by extact_asembly_stats.py',
                        type=str)
    io_opts.add_argument("--read-metadata",
                        dest="reads",
                        required=True,
                        help='Isolate read metadata output by extact_read_metadata.py',
                        type=str)
    args = parser.parse_args()
    return (args)

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()
    sys.stderr.write("\nImporting metadata files\n")
    with open(args.assemblies, "r") as metaFile:
        metadata = json.loads(metaFile.read())["information"]
    with open(args.reads, "r") as metaFile:
        metadata += json.loads(metaFile.read())["information"]
    sys.stderr.write("\nIterating through isolates and extracting assembly stats\n")
    data = {}
    for isolate in tqdm(metadata):
        if "coli" in isolate["Organism_name"]:
            isolate_stats = {}
            isolate_stats["contig_N50"] = isolate["contig_stats"]["N50"]
            isolate_stats["gc_content"] = isolate["contig_stats"]["gc_content"]
            isolate_stats["contig_mean_length"] = isolate["contig_stats"]["mean"]
            isolate_stats["contig_count"] = isolate["contig_stats"]["sequence_count"]
            isolate_stats["genome_length"] = isolate["contig_stats"]["total_bps"]
            isolate_stats["scaffold_N50"] = isolate["contig_stats"]["N50"]
            isolate_stats["scaffold_mean_length"] = isolate["scaffold_stats"]["mean"]
            isolate_stats["scaffold_count"] = isolate["scaffold_stats"]["sequence_count"]
            try:
                data[isolate["isolateNameUnderscore"]] = isolate_stats
            except KeyError:
                data[isolate["isolateName"]] = isolate_stats
    with open("ESC_population_assembly_stats.json", "w") as outJSON:
        outJSON.write(json.dumps(data))
    sys.exit(0)

if __name__ == '__main__':
    main()