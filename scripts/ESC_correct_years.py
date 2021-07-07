import json
from tqdm import tqdm

with open("retrieved_ena_read_metadata/isolateReadAttributes.json") as inJson:
    metadata = json.loads(inJson.read())["information"]

for line in tqdm(metadata):
    if not "Year" in line and "BioSample_CollectionDate" in line:
        try:
            year = line["BioSample_CollectionDate"].split("-")[0]
            line["Year"] = int(year)
        except ValueError:
            print(line["BioSample_CollectionDate"])

with open("retrieved_ena_read_metadata/isolateReadAttributes.json", "w") as outJson:
    outJson.write(json.dumps({"information": metadata}))