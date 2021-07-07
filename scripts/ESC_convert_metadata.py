from io import StringIO
import json
import pandas as pd
import sys
from tqdm import tqdm

sys.stderr.write("\nExtracting metadata from ESC csv file\n")
metadata = pd.read_csv("F1_genome_metadata.csv")

with open("FILTERED_MD_FINAL_ALL.tab", "rb") as accessFile:
    additions = accessFile.read().splitlines()[1:]

Ids = []
accessions = []
for line in additions:
    line = line.decode('windows-1252')
    Ids.append(line.split("\t")[0].upper())
    accessions.append(line.split("\t")[15])
ERS_dict = {}
ENA_accessions = []
unkown_accessions = []
for index, row in tqdm(metadata.iterrows()):
    if not row["ID"] == "":
        if "ESC_" in row["ID"]:
            try:
                accession = accessions[list(Ids).index(row["ID"])]
            except ValueError:
                accession = row["ID"]
        else:
            accession = row["ID"]
        ERS_dict[accession] = {"Assembly_name": row["Assembly_name"],
                                "PopPUNK": row["PopPUNK"],
                                "ST": row["ST"],
                                "MDR": row["MDR"],
                                "Ab_classes": row["Ab_classes"],
                                "Pathotype": row["Pathotype"],
                                "Phylogroup": row["Phylogroup"],
                                "Isolation": row["Isolation"],
                                "Country": row["Country"],
                                "Continent": row["Continent"],
                                "Total_AMR_genes": row["Total_AMR_genes"],
                                "Total_virulence_genes": row["Total_virulence_genes"],
                                "Pathotype_Vir_genes": row["Pathotype_Vir_genes"],
                                "Other_Vir_genes": row["Other_Vir_genes"],
                                "AMR_genes": row["AMR_genes"]}
        if "ESC_" in accession or "," in accession:
            unkown_accessions.append(accession)
        else:
            ENA_accessions.append(accession)

with open("ESC_METADATA.json", "w") as m:
    m.write(json.dumps(ERS_dict))
with open("ESC_ENA_Accessions.txt", "w") as e:
    e.write("\n".join(ENA_accessions))
with open("ESC_UNKNOWN_Accessions.txt", "w") as e:
    e.write("\n".join(unkown_accessions))

sys.stderr.write("Done\n")
