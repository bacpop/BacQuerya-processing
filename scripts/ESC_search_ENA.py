"""Search for term in EBI to get accessions that link to genomic metadata in XML format"""
from joblib import Parallel, delayed
import requests
from tqdm import tqdm
import xmltodict

with open("ENA_READS_SUPPRESSED.txt", "r") as inAccess:
    Ids = inAccess.read().splitlines()

def search_ENA(Id):
    if Id == "ND":
        return None
    apiURL = "https://www.ebi.ac.uk/ebisearch/ws/rest/sra-sample?query={}&fields=acc,description,name&size=100".format(Id.replace("#", "%23"))
    urlResponse = requests.get(apiURL)
    accession_metadata = dict(xmltodict.parse(urlResponse.text))
    try:
        if accession_metadata["result"]["hitCount"] == "0":
            apiURL = "https://www.ebi.ac.uk/ebisearch/ws/rest/wgs_masters?query={}&fields=acc,description,name&size=100".format(Id.replace("#", "%23"))
            urlResponse = requests.get(apiURL)
            accession_metadata = dict(xmltodict.parse(urlResponse.text))
            if accession_metadata["result"]["hitCount"] == "0":
                if "NC_" in Id:
                    with open("ESC_NCBI_ASSEMBLIES.txt", "a") as assem:
                        assem.write(Id + "\n")
                else:
                    with open("ESC_NOT_FOUND.txt", "a") as suppressed:
                        suppressed.write(Id + "\n")
            else:
                try:
                    access = accession_metadata["result"]["entries"]['entry']["@id"]
                except:
                    if "NC_" in Id:
                        with open("ESC_NCBI_ASSEMBLIES.txt", "a") as assem:
                            assem.write(Id + "\n")
                    else:
                        with open("ESC_NOT_FOUND.txt", "a") as suppressed:
                            suppressed.write(Id + "\n")
                txtURL = "https://www.ebi.ac.uk/ena/browser/api/embl/{}?lineLimit=1000".format(access)
                urlResponse = requests.get(txtURL).text.splitlines()
                for resp in urlResponse:
                    if "BioSample" in resp:
                        biosample = resp.split(" ")[4].replace(".", "")
                        with open("ESC_FOUND.txt", "a") as found:
                            found.write(Id + "\n")
                        with open("ESC_SEARCH_FOUND_ACCESSIONS.txt", "a") as identified:
                            identified.write(biosample + "\n")
                        return biosample
        else:
            access = accession_metadata["result"]["entries"]['entry']["@id"]
            with open("ESC_FOUND.txt", "a") as found:
                found.write(Id + "\n")
            with open("ESC_SEARCH_FOUND_ACCESSIONS.txt", "a") as identified:
                identified.write(access + "\n")
            return access
    except:
        with open("ESC_NOT_FOUND.txt", "a") as suppressed:
            suppressed.write(Id + "\n")

job_list = [
        Ids[i:i + 20] for i in range(0, len(Ids), 20)
    ]
foundAccessions = []
for job in tqdm(job_list):
    foundAccessions +=  Parallel(n_jobs=20)(delayed(search_ENA)(i) for i in job)