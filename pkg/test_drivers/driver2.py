import dtools.cluster as sl
import numpy as np

import pkg.dtools.misc

n = 200
L = 10
xx = pkg.dtools.misc.simulate(lambda x: x ** 2, 5, n, L)
yy = pkg.dtools.misc.simulate(lambda x: -x, 5, 200, 10)
zz = pkg.dtools.misc.simulate(lambda x: x, 5, 200, 10)
data = np.concatenate((xx, yy, zz))
colors = ['r'] * 200 + ['b'] * 200 + ['g'] * 200


def eval2(a, b, N):
    T = 0.0
    total = 0.0
    for i in xrange(N):
        for j in xrange(i+1, N):
           if (a[i] == a[j] and b[i] == b[j]) or\
                (a[i] != a[j] and b[i] != b[j]):
               T += 1
           total += 1
    return T/total


d = sl.Clusterer(data)
d.approximate(4,'dwt')
d.perform_sax('abc')
d.saxcluster(lookup='False')
