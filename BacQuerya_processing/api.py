import cobs_index as cobs
from flask import Flask, request, jsonify
from flask_cors import CORS, cross_origin
import h5py
import json
import os
import pp_sketchlib
import sys
from types import SimpleNamespace

sys.path.insert(1, '..')

# allow upload of entire assembly
# post entire assembly to backend and save
# create pp-sketchlib database
# query database
# return most closely related species name
# query the correct cobs index

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

def constructQueryDB(q_file, query_db):
    names = []
    sequences = []
    with open(q_file, 'r') as queryFile:
        for queryLine in queryFile:
            queryFields = queryLine.rstrip().split("\t")
            names.append(queryFields[0])
            sequences.append(list(queryFields[1:]))
    pp_sketchlib.constructDatabase(query_db, names, sequences, [21, 23],
                                   int(round(1000/64)),
                                   False,
                                   False,
                                   True,
                                   20,
                                   False,
                                   3,
                                   False,
                                   0)

def getSampleNames(db_prefix):
    rList = []
    ref = h5py.File(db_prefix + ".h5", 'r')
    for sample_name in list(ref['sketches'].keys()):
        rList.append(sample_name)
    return(rList)

def iterDistRows(refSeqs, querySeqs, self=True):
    """Gets the ref and query ID for each row of the distance matrix
    Returns an iterable with ref and query ID pairs by row.
    Args:
        refSeqs (list)
            List of reference sequence names.
        querySeqs (list)
            List of query sequence names.
        self (bool)
            Whether a self-comparison, used when constructing a database.
            Requires refSeqs == querySeqs
            Default is True
    Returns:
        ref, query (str, str)
            Iterable of tuples with ref and query names for each distMat row.
    """
    if self:
        if refSeqs != querySeqs:
            raise RuntimeError('refSeqs must equal querySeqs for db building (self = true)')
        for i, ref in enumerate(refSeqs):
            for j in range(i + 1, len(refSeqs)):
                yield(refSeqs[j], ref)
    else:
        for query in querySeqs:
            for ref in refSeqs:
                yield(ref, query)

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

@app.route('/assembly', methods=['POST'])
@cross_origin()
def analyseAssembly():
    if request.method == 'POST':
        query_dir = "uploaded_assembly"
        if not os.path.exists(query_dir):
            os.mkdir(query_dir)
        query_filename = os.path.join(query_dir, "q_file.txt")
        uploaded_file = request.files['file']
        uploaded_file.save(os.path.join(query_dir, 'query_file.txt'))
        with open(query_filename, "w") as w:
            w.write("query\t" + os.path.join(query_dir,"query_file.txt"))
        # construct pp-sketchlib query db
        query_db = os.path.join(query_dir, query_dir)
        constructQueryDB(query_filename, query_db)
        # query genome browser database to identify species
        ref_db = "genome_browser/taxons/taxons"
        rList = getSampleNames(ref_db)
        qList = getSampleNames(query_db)
        ref = h5py.File(ref_db + ".h5", 'r')
        query = h5py.File(query_db + ".h5", 'r')
        db_kmers = list(set(ref['sketches/' + rList[0]].attrs['kmers']).intersection(
           query['sketches/' + qList[0]].attrs['kmers']
        ))
        distMat = pp_sketchlib.queryDatabase(ref_db, query_db, rList, qList, db_kmers,
                                             True, False, 1, False,
                                             0)
        names = iterDistRows(rList, qList, rList == qList)
        sys.stdout.write("\t".join(['Query', 'Reference', 'Core', 'Accessory']) + "\n")
        for i, (ref, query) in enumerate(names):
            sys.stdout.write("\t".join([query, ref, str(distMat[i,0]), str(distMat[i,1])]) + "\n")
        return jsonify({"test": "response"})

if __name__ == "__main__":
    app.secret_key = os.urandom(24)
    app.run(debug=False,use_reloader=False)

