import cobs_index as cobs
from flask import Flask, request, jsonify
from flask_cors import CORS, cross_origin
import os
import sys

sys.path.insert(1, '..')

app = Flask(__name__)
CORS(app, expose_headers='Authorization')

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
def postJSON():
    if not request.json:
        return "not a json post"
    if request.json:
        sequence_dict = request.json
        query_sequence = sequence_dict['searchTerm']

        index = cobs.Search()
        result = index.search(sequence, threshold = 0.8)
        response = {"sequenceResult" : result}
    return jsonify(response)

if __name__ == "__main__":
    app.secret_key = os.urandom(24)
    app.run(debug=False,use_reloader=False)

