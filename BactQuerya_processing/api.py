import cobs_index as cobs
from flask import Flask, request, jsonify
from flask_cors import CORS, cross_origin
import json
import os
import sys

sys.path.insert(1, '..')

app = Flask(__name__)
CORS(app, expose_headers='Authorization')

def load_metrics(match_count,
                 index,
                 gene_dicts_list,
                 query_length,
                 kmer_length):
    """Search through metadata to find corresponding index"""
    for gene_dict in gene_dicts_list:
        isolate = gene_dict["isolateName"]
        feature_list = gene_dict["features"]
        for feature in feature_list:
            if "gene_index" in feature.keys():
                if int(feature["gene_index"]) == index:
                    match_proportion = round(match_count*100/((query_length-kmer_length)+1), 2)
                    result_metrics = {"isolateName": isolate,
                                    "geneName": feature["Name"][0],
                                    "numberMatching": match_proportion,
                                    "gene_index": index}
    return result_metrics

@app.route('/JSON', methods=['POST'])
@cross_origin()
def postJSON():
    if not request.json:
        return "not a json post"
    if request.json:
        with open("test/test_isolates/isolateFeatures.json", "r") as j:
            json_file = j.read()
    return jsonify(json_file)

@app.route('/sequence', methods=['POST'])
@cross_origin()
def postSeqResult():
    if not request.json:
        return "not a json post"
    if request.json:
        sequence_dict = request.json
        query_sequence = sequence_dict['searchTerm']
        # search for uploaded sequence in COBS index
        index_name = "sequence_index/20_index.cobs_compact"
        index = cobs.Search(index_name)
        result = index.search(query_sequence, threshold = 0.8)
        # load metadata for identified sequences
        with open("isolate_genes/allIsolates.json") as f:
            geneJSON = f.read()
        gene_dicts_list = json.loads(geneJSON)
        gene_dicts_list = gene_dicts_list["information"]
        query_length = len(query_sequence)
        kmer_length = int(os.path.basename(index_name).split("_")[0])
        result_metrics = []
        for res in result:
            match_count = res[0]
            index = int(res[1])
            metrics = load_metrics(match_count,
                                   index,
                                   gene_dicts_list,
                                   query_length,
                                   kmer_length)
            result_metrics.append(metrics)
        response = {"resultMetrics" : result_metrics}
    return jsonify(response)

if __name__ == "__main__":
    app.secret_key = os.urandom(24)
    app.run(debug=False,use_reloader=False)

