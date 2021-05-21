"""Convert GPS metadata csv files to a json for later input into extract_assembly_stats.py and extract_read_metadata.py"""
import json
import os
import pandas as pd
import sys
from tqdm import tqdm

sys.stderr.write("\nExtracting metadata from GPS csv files\n")
# if this is the GPS data, we are going to supplement the metadata with those in the mmc supplementary data and GPS1_in_silico_output.csv
GPS_insilicoDF = pd.read_csv("GPS1_insilico_output.csv")
GPS_insilicoDict = {}
for index, row in GPS_insilicoDF.iterrows():
    if not row["Lane_Id"] == "":
        GPS_insilicoDict[row["Lane_Id"]] = {"ERS": row["ERS"],
                                            "ERR": row["ERR"],
                                            "Lane_Id": row["Lane_Id"],
                                            "In_Silico_Serotype": row["in_siico_serotype"],
                                            "In_Silico_St": row["In_Silico_St"],
                                            "WGS_PEN_SIR_Nonmeningitis": row["WGS_PEN_SIR_Nonmeningitis"]}
del GPS_insilicoDF
GPS_mmc1 = pd.read_csv("mmc_tab1.csv")
GPS_supplement = {}
for index, row in tqdm(GPS_mmc1.iterrows()):
    supplement_dict = {"ERR": row["ERR"],
                        "Lane_Id": row["ID"],
                        "In_Silico_St": row["In_Silico_St"],
                        "CC": row[" CC"],
                        "In_Silico_Serotype": row["In_Silico_Serotype"],
                        "GPSC": row["GPSC"],
                        "GPSC_PCV_type_pre_PCV": row["GPSC-PCV-type pre-PCV"],
                        "Vaccine_Period": row["Vaccine_Period"],
                        "Vaccine_Status": row["Vaccine_Status"],
                        "Country": row["Country"],
                        "Year": row["Year"],
                        "Date(month)": row["Date (mid-year if day/month unknown)"],
                        "Age_group": row["Age_group"],
                        "Disease": row["Disease"]}
    if row["ID"] in GPS_insilicoDict.keys():
        supplement_dict.update(GPS_insilicoDict[row["ID"]])
    GPS_supplement[supplement_dict["ERR"]] = supplement_dict
# if not in mmc1 then add manually
for key, value in GPS_insilicoDict.items():
    if not value["ERR"] in GPS_supplement.keys():
        GPS_supplement[value["ERR"]] = value
with open("GPS_metadata.json", "w") as outMeta:
    outMeta.write(json.dumps(GPS_supplement))
