import subprocess
import sys

# download entries
sys.stderr.write("\nRetrieving genomic sequences (-a genome)\n")
subprocess.run("python ../entrez_extract-runner.py -s test_accessions.txt -e danielanderson1@hotmail.com --threads 3 -o test_files -a genome", shell=True, check=True)

sys.stderr.write("\nRetrieving annotations (-a annotation)\n")
subprocess.run("python ../entrez_extract-runner.py -s test_accessions.txt -e danielanderson1@hotmail.com --threads 3 -o test_files -a annotation", shell=True, check=True)

sys.stderr.write("\nRetrieving assembly stats (-a assembly-stats)\n")
subprocess.run("python ../entrez_extract-runner.py -s test_accessions.txt -e danielanderson1@hotmail.com --threads 3 -o test_files -a assembly-stats", shell=True, check=True)

sys.stderr.write("\nRetrieving assembly report (-a assembly-report)\n")
subprocess.run("python ../entrez_extract-runner.py -s test_accessions.txt -e danielanderson1@hotmail.com --threads 3 -o test_files -a assembly-report", shell=True, check=True)

# unqip genome and annotation
subprocess.run("gunzip test_files/*.gz", shell=True, check=True)

# extract features from entries
sys.stderr.write("\nCoverting GFF files to JSON strings\n")
subprocess.run("python ../feature_extract-runner.py -s test_files -g test_files -o test_features --threads 3", shell=True, check=True)
