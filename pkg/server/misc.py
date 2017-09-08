"""
This stores the global variables
"""
import cPickle as pickle
import collections
import random

cache = collections.defaultdict(dict)

def get_data(user, token, path):
    global cache
    if user not in cache or token not in cache[user] or not token:
        token = "%032x" % random.getrandbits(128)
        with open(path, 'rb') as f:
            cache[user][token] = pickle.load(f)
            f.close()
    return cache[user][token],  token


def get_obj(user, token, path):
    global obj_cache
    if user not in cache or token not in cache[user] or not token:
        token = "%032x" % random.getrandbits(128)
        with open(path, 'rb') as f:
            cache[user][token] = pickle.load(f)
            f.close()
    return cache[user][token],  token


def clear_cache(user):
    global cache
    if user in cache:
        del cache[user]
    if user in obj_cache:
        del obj_cache[user]


def convert(input):
    if isinstance(input, dict):
        return {convert(key): convert(value) for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [convert(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input