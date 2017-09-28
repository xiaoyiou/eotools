import collections
import hashlib
import json

import bencode
from flask import Flask, request, make_response, abort, jsonify

import misc
from pkg.core import mdata
from ..core import rawdata

app = Flask(__name__)

obj_cache = misc.obj_cache
# some helper functions
def get_field(x, y, req): return req[x] if x in req else y


def get_fields(xs, ys, req):
    assert len(xs) == len(ys)
    return [get_field(x, y, req) for x, y in zip(xs, ys)]



@app.route('/')
def index():
    return "Hello, World!"

@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify({'error': 'Not found'}), 404)


@app.route('/api/cexp', methods=['POST'])
def create_exp():
    if app.debug:
        print request.json
    # required fields in request
    req = json.loads(request.json)
    fields = ['tpaths', 'names', 'clpath', 'gpath', 'path']
    if not req or any([not x in req for x in fields]):
        abort(400)
    if len(req['tpaths']) != len(req['names']):
        abort(400)

    # retrieving fields
    keys = fields + ['us', 'ds', 'overlap', 'method', 'type', 'power', 'window']
    defaults = [''] * 5 + [1000, 1000, 0.5, 'mean', 'bed', 4, 100]
    paths, names, clpath, gpath, path, us, ds, overlap, method, type, power, window = get_fields(keys, defaults, req)
    exp = rawdata.Experiment(paths, names)

    # creating experiment object
    exp.preprocess(clpath, gpath, window, us, ds, overlap, method, 4, type, n_workers=8)
    exp.get_binary(power=power)
    exp.export(req['path'])

    return make_response(jsonify({'result': 'completed'}), 200)


@app.route('/api/getraw', methods=['POST'])
def get_raw():
    if app.debug:
        print request.json
    req = json.loads(request.json)
    fields = ['user', 'glst', 'mlst', 'path']
    if not req or any([not x in req for x in fields]):
        abort(400)
    keys = fields + ['token']
    defaults = [''] * 4 + [None]
    user, glst, mlst, path, token = get_fields(keys, defaults, req)

    # trying to get cache
    exp, token = misc.get_data(user, token, path)

    ret = {}
    genes = set(exp.genes)
    for mod in mlst:
        tmp = {}
        subset = exp.raw[mod].loc[list(genes.intersection(glst))]
        for ind, row in subset.iterrows():
            tmp[ind] = row.tolist()
        ret[mod] = tmp
    ret['token'] = token
    return make_response(json.dumps(ret, ensure_ascii=False, indent=4), 200)


@app.route('/api/getpattern', methods=['POST'])
def get_pattern():

    if app.debug:
        print request.json
    req = json.loads(request.json)
#    req = request.json
    global obj_cache
    fields = ['user', 'label']
    if not req or any([x not in req for x in fields]):
        abort(400)

    keys = fields + ['token', 'labels', 'topn', 'path']
    defaults = [''] * 2 + [None, None, 100, '']
    user, label, token, labels, topn, path = get_fields(keys, defaults, req)
    if token is None:
        if labels is None:
            abort(400)
        token = hashlib.md5(bencode.bencode([misc.convert(path),
                                             misc.convert(labels)])).hexdigest()
    aobj = None
    if user in obj_cache and token in obj_cache[user]:
        aobj = obj_cache[user][token]
    else:
        exp, _ = misc.get_data(user, None, path)
        aobj = mdata.Mdata(exp.binary, exp.names)
        aobj.addLabels(exp.genes, 'all')
        genes = set(exp.genes)
        for label in labels:
            lst = list(genes.intersection(labels[label]))
            if len(lst) == 0:
                abort(400)
            aobj.addLabels(lst, label)
        
        aobj.findPatterns(2)
        aobj.findPScores(mode='ar', ratios=aobj.ratios)
        obj_cache[user][token] = aobj
    ret = {'token': token, 'patterns': []}
    for ind, row in aobj.pScore.sort_values(by=[label], ascending=False).head(topn).iterrows():
        ret['patterns'].append({"pattern": list(ind), "score":row[label]})
        ret['tnames'] = aobj.mods
    return make_response(json.dumps(ret, ensure_ascii=False, indent=4), 200)

@app.route('/api/getpgenes', methods=['POST'])
def get_pgenes():
    if app.debug:
        print request.json
    req = json.loads(request.json)
    fields = ['user', 'token', 'pattern', 'label']
    if not req or any([not x in req for x in fields]):
        abort(400)
    defualts = [None] * 4
    user, token, pattern, label = get_fields(fields, defualts, req)
    if user not in obj_cache or token not in obj_cache[user]:
        abort(400)
    pattern = frozenset(map(int, get_field('pattern', None, req)))
    aobj = obj_cache[user][token]
    filt = set(aobj.labels[aobj.labels[label] == 1].index.tolist())
    lst = list(filt.intersection(aobj.pgMap[aobj.pgMap[pattern] == 1].index.tolist()))
    
    return make_response(json.dumps({"genes": lst}, ensure_ascii=False, indent=4), 200)



@app.route('/api/clear_cache', methods=['POST'])
def clear_cache():
    if app.debug:
        print request.json
    req = json.loads(request.json)
    fields = ['user']
    if not req or any([not x in req for x in fields]):
        abort(400)
    user = req['user']
    global obj_cache
    if user in obj_cache:
        del obj_cache[user]
    misc.clear_cache(user)
    return make_response(jsonify({'result': 'completed'}), 200)
    
