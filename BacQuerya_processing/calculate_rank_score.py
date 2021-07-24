"""A script to calculate a ranking score for isolates based on the available metadata. The score is used as a tie break when querying the BacQuerya elastic indices."""
from joblib import Parallel, delayed
import json
import numpy as np
import sys
from tqdm import tqdm

def get_options():

    import argparse

    parser = argparse.ArgumentParser(prog='calculate_rank_score')
    io_opts = parser.add_argument_group('input')
    io_opts.add_argument("--assembly-metadata",
                        dest="assemblies",
                        required=False,
                        default=None,
                        help='Isolate asssembly metadata output by extact_asembly_stats.py',
                        type=str)
    io_opts.add_argument("--read-metadata",
                        dest="reads",
                        required=False,
                        default=None,
                        help='Isolate read metadata output by retrieve_read_metadata.py',
                        type=str)
    io_opts.add_argument("--threads",
                        dest="n_cpu",
                        required=False,
                        help="number of threads for extracting features",
                        default=1,
                        type=int)
    args = parser.parse_args()
    return (args)

def calculate_sdi(hashDict):
    """Calculate simpsons diversity index of the sample"""
    hashCounts = hashDict.values()
    hashCounts = [int(x.split("/")[0]) for x in hashCounts]
    total = sum(hashCounts)
    sdi = 0
    for count in hashDict.values():
        abundance = int(count.split("/")[0])/total
        sdi += abundance ** 2
    return 1 - sdi

def assess_contamination(hashDict, organism):
    """Determine if a genomic sequence is contaminated and adjust rank score. If SHI is above 0.2, -10 from rank score."""
    sdi = calculate_sdi(hashDict)
    contaminated = False
    for species, hashes in hashDict.items():
        if "seqs]" in species:
            if int(hashes.split("/")[0]) > 10 and not species.split(" ")[2].lower() in organism:
                contaminated = True
        else:
            if int(hashes.split("/")[0]) > 10 and not species.split(" ")[1].lower() in organism:
                contaminated = True
    if sdi > 0.2 or contaminated:
        return -10
    else:
        return 0

def caculate_score(isolate, GC_lower, GC_upper, length_lower, length_upper):
    """Calculate a rank score for all isolates"""
    score = 0
    # currently it takes too long to download read sets so we are only running mash on assemblies
    if "mashHashes" in isolate.keys():
        # calculate contamination and adjust score
        mashHashes = isolate["mashHashes"]
        mashSpecies = isolate["mashSpecies"]
        organism = isolate["Organism_name"]
        hashDict = {}
        for species in range(len(mashSpecies)):
            hashDict[mashSpecies[species]] = mashHashes[species]
        score += assess_contamination(hashDict, organism)
    # adjust score depending on if reads or assembly
    if isolate["Genome_representation"] == "reads":
        score += 1
    if isolate["Genome_representation"] == "full":
        score += 10
        # add 100 if 1 contig
        if int(isolate["contig_stats"]["sequence_count"]) == 1:
            score += 100
        # adjust score if assembly statistics show the isolate is an outlier
        if float(isolate["contig_stats"]["gc_content"]) > GC_upper or float(isolate["contig_stats"]["gc_content"]) < GC_lower:
            score -= 5
        if int(isolate["contig_stats"]["total_bps"]) > length_upper or int(isolate["contig_stats"]["total_bps"]) < length_lower:
            score -= 5
        if int(isolate["contig_stats"]["longest"]) < 10000:
            score -= 5
        if int(isolate["contig_stats"]["N50"]) < 5000:
            score -= 5
    # add 1 point for every field in the isolate metadata
    score += len(isolate.values())
    isolate["rankScore"] = score
    return isolate

def main():
    """Main function. Parses command line args and calls functions."""
    args = get_options()
    sys.stderr.write("\nImporting metadata files\n")
    metadata = []
    if args.assemblies:
        with open(args.assemblies, "r") as assemblyMeta:
            assemblies = json.loads(assemblyMeta.read())["information"]
        metadata += [assem for assem in assemblies if assem is not None]
    if args.reads:
        with open(args.reads, "r") as readMeta:
            reads = json.loads(readMeta.read())["information"]
        metadata += [r for r in reads if r is not None]
    # calculate population wide IQR for GC content and length
    GC_proportions  = []
    lengths = []
    sys.stderr.write("\nCalculating assembly statistic ranges\n")
    for isolate_row in tqdm(metadata):
        if isolate_row["Genome_representation"] == "full":
            GC_proportions.append(isolate_row["contig_stats"]["gc_content"])
            lengths.append(isolate_row["contig_stats"]["total_bps"])
    GC75, GC25 = np.percentile(GC_proportions, [75 ,25])
    GC_lower = GC25 - 1.5*(GC75 - GC25)
    GC_upper = GC75 + 1.5*(GC75 - GC25)
    length75, length25 = np.percentile(GC_proportions, [75 ,25])
    length_lower = length25 - 1.5*(length75 - length25)
    length_upper = length75 + 1.5*(length75 - length25)
    # parallelise calculation of the rank score for assemblies
    if args.assemblies:
        job_list = [assemblies[i:i + args.n_cpu] for i in range(0, len(assemblies), args.n_cpu)]
        assembliesScored = []
        sys.stderr.write("\nCalculating assembly scores\n")
        for job in tqdm(job_list):
            assembliesScored += Parallel(n_jobs=args.n_cpu)(delayed(caculate_score)(metaLine,
                                                                                    GC_lower,
                                                                                    GC_upper,
                                                                                    length_lower,
                                                                                    length_upper) for metaLine in job)
        with open(args.assemblies, "w") as assembliesOut:
            assembliesOut.write(json.dumps({"information": assembliesScored}))
    # parallelise calculation of the rank score for reads
    if args.reads:
        job_list = [reads[i:i + args.n_cpu] for i in range(0, len(reads), args.n_cpu)]
        readsScored = []
        sys.stderr.write("\nCalculating read scores\n")
        for job in tqdm(job_list):
            readsScored += Parallel(n_jobs=args.n_cpu)(delayed(caculate_score)(metaLine,
                                                                            GC_lower,
                                                                            GC_upper,
                                                                            length_lower,
                                                                            length_upper) for metaLine in job)
        with open(args.reads, "w") as readsOut:
            readsOut.write(json.dumps({"information": readsScored}))
    sys.exit(0)

if __name__ == '__main__':
    main()