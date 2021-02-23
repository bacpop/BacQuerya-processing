import subprocess
import sys


def get_options():
    import argparse
    parser = argparse.ArgumentParser(prog='run_test')
    e = parser.add_argument_group('email')
    e.add_argument("-e",
                    "--email",
                    dest="email",
                    required=True,
                    help='user email address for entrez access',
                    type=str)
    args = parser.parse_args()
    return (args)

args = get_options()

# download entries
sys.stderr.write("\nRetrieving genomic sequences (-a genome)\n")
subprocess.run("python ../entrez_extract-runner.py -s test_accessions.txt -e " + args.email + " --threads 3 -o test_files -a genome", shell=True, check=True)

sys.stderr.write("\nRetrieving annotations (-a annotation)\n")
subprocess.run("python ../entrez_extract-runner.py -s test_accessions.txt -e " + args.email + " --threads 3 -o test_files -a annotation", shell=True, check=True)

sys.stderr.write("\nRetrieving assembly stats (-a assembly-stats)\n")
subprocess.run("python ../entrez_extract-runner.py -s test_accessions.txt -e " + args.email + " --threads 3 -o test_files -a assembly-stats", shell=True, check=True)

sys.stderr.write("\nRetrieving assembly report (-a assembly-report)\n")
subprocess.run("python ../entrez_extract-runner.py -s test_accessions.txt -e " + args.email + " --threads 3 -o test_files -a assembly-report", shell=True, check=True)

# extract attributes from assembly stats
sys.stderr.write("\nCoverting Assembly stats to JSON strings\n")
subprocess.run("python ../isolate_attributes-runner.py -a test_files -i 0 -o test_attributes/isolateFeatures.json --threads 3", shell=True, check=True)

# extract features from entries
sys.stderr.write("\nCoverting GFF files to JSON strings\n")
subprocess.run("python ../feature_extract-runner.py -s test_files -g test_files -j test_attributes/isolateFeatures.json -o test_features/allGenes.json --threads 3", shell=True, check=True)
