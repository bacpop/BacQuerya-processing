import cobs_index as cobs
from flask import Flask, request, jsonify
from flask_cors import CORS, cross_origin
import h5py
import json
import os
import pp_sketchlib
import re
import sys
from tqdm import tqdm
from types import SimpleNamespace

from paper_search import search_pubmed

sys.path.insert(1, '..')

# allow upload of entire assembly
# post entire assembly to backend and save
# create pp-sketchlib database
# query database
# return most closely related species name
# query the correct cobs index

app = Flask(__name__)
CORS(app, expose_headers='Authorization')

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
        sys.stderr.write("\nSearching COBS index\n")
        index_name = "index_genes/31_index.cobs_compact"
        index = cobs.Search(index_name)
        result = index.search(query_sequence, threshold = 0.8)
        # load metadata for identified sequences
        sys.stderr.write("\nLoading gene metadata\n")
        with open("extracted_genes/panarooPairs.json") as f:
            geneJSON = f.read()
        genePairs = json.loads(geneJSON)
        query_length = len(query_sequence)
        kmer_length = int(os.path.basename(index_name).split("_")[0])
        result_metrics = []
        sys.stderr.write("\nExtracting metadata for COBS result\n")
        for res in tqdm(result):
            match_count = int(res.score)
            index = res.doc_name.split("_")[0]
            for k, v in genePairs.items():
                if v == int(index):
                    geneName = k
            match_proportion = round(match_count*100/((query_length-kmer_length)+1), 2)
            metrics = {"geneName": geneName, "numberMatching": match_proportion}
            result_metrics.append(metrics)
        sys.stderr.write("\nPosting results to frontend\n")
        response = {"resultMetrics" : result_metrics}
    return jsonify(response)

@app.route('/species', methods=['POST'])
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

@app.route('/assembly', methods=['GET, POST'])
@cross_origin()
def searchAssemblyIndex():
    if request.method == 'POST':
        print("File uploaded!")
        query_dir = "uploaded_assembly"
        if not os.path.exists(query_dir):
            os.mkdir(query_dir)
        uploaded_file = request.files['file']
        uploaded_file.save(os.path.join(query_dir, 'query_file.txt'))
        with open(os.path.join(query_dir, 'query_file.txt'), "r") as f:
            query_sequence = f.read()
        # search for uploaded assembly in COBS index
        index_name = "index_assemblies/31_index.cobs_compact"
        index = cobs.Search(index_name)
        query_sequence = re.sub(r'[^ACTG]', '', query_sequence)
        query_sequence = query_sequence[:3000000]
        result = index.search(query_sequence, num_results = 10)
        # calculate match proportion
        query_length = len(query_sequence)
        resultList = []
        for match in result:
            match_count = int(match.score)
            matchingIsolate = match.doc_name
            kmer_length = int(os.path.basename(index_name).split("_")[0])
            match_proportion = round(match_count*100/((query_length-kmer_length)+1), 2)
            resultList.append({"isolateName": matchingIsolate,
                               "matchCount": match_count,
                               "matchProportion": match_proportion})
        return jsonify({"result" : resultList})
    else:
        print("No file uploaded!")

from urllib.parse import unquote

@app.route('/paper', methods=['POST'])
@cross_origin()
def paperSearch():
    if not request.json:
        return "not a json post"
    if request.json:
        searchDict = request.json
        searchTerm = searchDict["searchTerm"]
        maxResults = 100
        if searchDict["source"] == "URL":
            maxResults = 1
            searchTerm = unquote(searchTerm)
        searchResult = search_pubmed(searchTerm,
                                     "",
                                     maxResults)
        return jsonify({"result": searchResult})

if __name__ == "__main__":
    app.secret_key = os.urandom(24)
    app.run(debug=False,use_reloader=False)

