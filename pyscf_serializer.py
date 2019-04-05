import pickle 
import sys
from pyscf import gto, scf
import tempfile
import copy

class SerializerStream(object):
    pass

def eliminate_streams(obj, inplace=False):
    if isinstance(obj, file) or isinstance(obj, tempfile._TemporaryFileWrapper):
        return SerializerStream()
    if not inplace:
        obj = copy.deepcopy(obj)
    try:
        for attr in obj.__dict__.keys():
            if isinstance(getattr(obj, attr), file):
                setattr(obj, attr, SerializerStream())
            elif isinstance(getattr(obj, attr), tempfile._TemporaryFileWrapper):
                getattr(obj, attr).close()
                setattr(obj, attr, SerializerStream())
            else:
                setattr(obj, attr, eliminate_streams(getattr(obj, attr)))
    except AttributeError:
        pass
    return obj

def replace_streams(obj):
    if isinstance(obj, SerializerStream):
        return sys.stdout
    try:
        for attr in obj.__dict__.keys():
            if attr=='stdout':
                setattr(obj, attr, sys.stdout)
            elif hasattr(getattr(obj, attr), '__dict__'):
                replace_streams(getattr(obj, attr))
    except AttributeError:
        pass
    return obj

def serialize(obj, pfile='pyscf.pkl'):
    with open(pfile, 'wb') as f:
        pickle.dump(obj, f)

def deserialize(pfile='pyscf.pkl'):
    with open(pfile, 'rb') as f:
        return pickle.load(f)

def selective_serialize(obj, paths, pfile='pyscf.pkl', inplace=False):
    data = {}
    for path in paths:
        if isinstance(path, str):
            path = tuple(path.split('.'))
        tip = obj
        try:
            for node in path:
                tip = getattr(tip, node)
            data[path] = eliminate_streams(tip, inplace)
        except AttributeError:
            print 'Attribute "{}" not found.'.format(node)
            pass
    serialize(data, pfile)


def selective_deserialize(obj, pfile='pyscf.pkl'):
    data = deserialize(pfile)
    for path, datum in data.iteritems():
        tip = obj
        for node in path[:-1]:
            tip = getattr(tip, node)
        setattr(tip, path[-1], replace_streams(datum))
    return obj


def print_object_tree(key, obj, indent_level=0, keys=[]):
    if not hasattr(obj, '__dict__'):
        print ('\t'*indent_level)+'{}: {}'.format(keys, obj)
        return
    else:
        print keys
    for k, v in obj.__dict__.iteritems():
        print_object_tree(k, v, indent_level+1, keys+[k])

def object_type_tree(obj, path=[]):
    pairs = []
    pairs.append((path, type(obj)))
    if hasattr(obj, '__dict__'):
        for k, v in obj.__dict__.iteritems():
            pairs.extend(object_type_tree(v, path+[k]))
    return pairs









