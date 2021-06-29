"""Sample one isolate per GPSC"""
import json
from tqdm import tqdm

with open("GPS_poppunk.txt", "r") as GPSFile:
    isolateList = GPSFile.read().splitlines()
with open("GPS_metadata.json", "r") as metafile:
    GPSmetadata = json.loads(metafile.read())

metadataIsolates = []
for accession, data in GPSmetadata.items():
    metadataIsolates.append(data["Lane_Id"])

GPSC_counts = {}
subsamples = []
total_count = 0
for row in isolateList:
    splitRow = row.split("\t")
    laneID = splitRow[0].split("_")
    laneID = "#".join(["_".join(laneID[0:2])] + [laneID[2]])
    total_count += 1
    if not int(splitRow[1]) in GPSC_counts.keys():
        GPSC_counts[int(splitRow[1])] = 1
        if laneID in metadataIsolates:
            subsamples.append(laneID)
    else:
        GPSC_counts[int(splitRow[1])] += 1

with open("GPS_samples.txt", "w") as sampleOut:
    sampleOut.write("\n".join(subsamples))