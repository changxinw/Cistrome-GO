import time
import json
from pkg_resources import resource_filename


def Info(infoStr):
    print "[%s] %s" % (time.strftime('%H:%M:%S'), infoStr)


def readJson(file_name):
    res = json.load(open(resource_filename('cistromego', 'data/%s' % file_name)))
    return res


def readCoordinate(assembly):
    genome = resource_filename('cistromego','data/%s_nodup.refseq' % assembly)
    coordinate = []
    with open(genome, "rU") as gfile:
        for line in gfile.readlines():
            coordinate.append((line.rstrip().split("\t")[1], int(line.rstrip().split("\t")[6]), 1, line.rstrip().split("\t")[0]))
    return coordinate


def calculateRank(vector, rev=False):
    a={}
    rank=1
    for num in sorted(vector, reverse=rev):
        if num not in a:
            a[num]=rank
            rank=rank+1
    return [a[i] for i in vector]


def calculateRank2(vector, rev=False):
    vector_indice_tuple = [(vector[i], i) for i in range(len(vector))]
    vector_indice_tuple.sort(key=lambda x: x[0], reverse=rev)
    indices = [0] * len(vector)
    for i, j in enumerate(vector_indice_tuple):
        indices[j[1]] = i+1
    return indices


def normMinMax(vector):
    c_min = min(vector)
    c_max = max(vector)
    c_gap = c_max - c_min
    res_vec = map(lambda x: (x-c_min)/c_gap, vector)
    return res_vec
