# this file is just used to test the module
import random
import time

import numpy as np

import pkg.dtools.cluster as sl

np.random.seed(10)
N, D = 2000, 10
now = time.time()
test = np.random.random((N, D))
a = sl.Clusterer(test)
a.dtwcluster()
print (time.time() - now)


def simulate(func, w,  n):
    means = [func(x) for x in range(-n,n+1)]
    xx = [x + random.randint(-w,w) for x in range(-n, n + 1)]
    y = np.asarray([func(x) for x in xx]) + np.random.normal(means, 1, 2*n + 1)
    return y.tolist()


def eval(a, b, N):
    T = 0.0
    for i in xrange(N):
        for j in xrange(i+1, N):
           if (a[i] == a[j] and b[i] == b[j]) or\
                (a[i] != a[j] and b[i] != b[j]):
               T += 1
    return T/N/(N-1)

times = []
perfs = []
lens = [2, 4, 8]
for l in lens:
    ts = []
    perf = []
    for alph in xrange(2,6):
        now = time.time()
        alphabet = map(str,range(alph))
        b = sl.Clusterer(test)
        b.approximate(l, 'dwt', 'haar')
        b.perform_sax(alphabet)
        b.saxcluster(lookup=False)
        ts.append(time.time() - now)
        perf.append(eval(a, b, N))
    times.append(ts)
    perfs.append(perf)