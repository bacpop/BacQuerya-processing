"""Add/ modify fields for documents in the Elastic isolate index by index number"""
from elasticsearch import Elasticsearch
from joblib import Parallel, delayed
import sys
from tqdm import tqdm

sys.path.append(".")
from BacQuerya_processing.secrets import ELASTIC_API_URL, ELASTIC_ISOLATE_API_ID, ELASTIC_ISOLATE_API_KEY, SQL_CONNECTION_STRING
from BacQuerya_processing.extract_assembly_stats import standardise_species

def modify_entry(index_no, supplementary_dict):
    """Update an elastic document from a dictionary"""
    client = Elasticsearch([ELASTIC_API_URL],
                            api_key=(ELASTIC_ISOLATE_API_ID, ELASTIC_ISOLATE_API_KEY))
    fetchData = {
            "query": {
                "terms": {
                "_id": [str(index_no)]
                }
            }
        }
    isolateResult = client.search(index = "isolate_index",
                                  body = fetchData,
                                  request_timeout = 20)
    if not len(isolateResult["hits"]["hits"]) == 0:
        result = isolateResult["hits"]["hits"][0]["_source"]
        #result.update(supplementary_dict)
        if "GenBank_assembly_accession" in result:
            result["GenBank_assembly_accession"] = result["GenBank_assembly_accession"].replace(" ", "_")
            response = client.index(index = "isolate_index",
                                    id = result["isolate_index"],
                                    body = result,
                                    request_timeout=60)

indices = [i for i in range(36992)]
job_list = [indices[i:i + 8] for i in range(0, len(indices), 8)]
for job in tqdm(job_list):
    Parallel(n_jobs=8)(delayed(modify_entry)(i, {}) for i in job)