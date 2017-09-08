"""
ADD SOMETHING
This module provides the clustering methods
and stores the results in the object.
Now it's using affinity propagation, maybe in the future we
can use other clustering methods.
"""
import sdata as st
import collections
import pandas as pd
from sklearn.cluster import AffinityPropagation
import numpy as np
from fastdtw import dtw as dtwf
from dtw import dtw as dtws
from scipy.spatial.distance import euclidean

class Clusterer(object):
    def __init__(self, data):
        self.data = data
        self.N = data.shape[0]
        self.alphabets = None
        self.n_parts = None
        self.lookup = None   # the distance matrix for SAX method
        self.indexes = None  #A the SAX groups
        self.dists = None
        self.bins = None
        self.avdata = dict()
        self.saxsize = dict()
        # The info about the clustering result
        self.n_clusters = 0
        self.clusters = None
        self.cluster_centers = None
        self.cluster_sax = None
        self.asax_data = None
        self.ass = None     # the cluster assignments for each data point, a long list ofa
                            # sax ids

    def approximate(self, n, mode=None, wavelet='db1', norm=True):
        if mode is None:
            return
        tmp = st.approximate(self.data, n, mode, wavelet)
        self.data = st.znorm(tmp) if norm else tmp
        self.n_parts = self.data.shape[1]

    def perform_sax(self, alphabets):
        self.alphabets = alphabets
        sax_data, self.lookup, self.bins = st.sax_transform(self.data, self.alphabets)
        temp = self.data
        self.indexes = pd.DataFrame(sax_data).groupby(by=range(temp.shape[1])).groups
        for sax in self.indexes:
            self.saxsize[sax] = len(self.indexes[sax])
            self.avdata[sax] = temp[self.indexes[sax], :].mean(axis=0)

    def saxcluster(self, preference=None, lookup=True):

        cls = AffinityPropagation(preference=preference, affinity='precomputed') if lookup else \
            AffinityPropagation(preference=preference)
        if self.dists is None:
            if lookup:
                data = self.dists = self.__saxDists()
            else:
                data = self.dists = self.avdata.values()
        else:
            data = self.dists
        cls.fit(data)
        reps = self.indexes.keys()
        self.cluster_sax = [reps[i] for i in cls.cluster_centers_indices_]
        self.cluster_centers = [self.avdata[sax] for sax in self.cluster_sax]
        self.clusters = collections.defaultdict(list)
        for ind, label in enumerate(cls.labels_):
            sax = self.cluster_sax[label]
            self.clusters[sax] += self.indexes.values()[ind]
        self.asax_data = dict()
        for sax in self.clusters:
            self.asax_data[sax] = self.data[self.clusters[sax], :].mean(axis=0)
        self.ass = [0] * self.N
        for sax in self.cluster_sax:
            v = self.cluster_sax.index(sax)
            for ind in self.clusters[sax]:
                self.ass[ind] = v
        self.n_clusters = len(self.clusters)

    def dtwcluster(self, window=None, fast=True):
        cls = AffinityPropagation(preference=None, affinity='precomputed')
        window = self.data.shape[1] if fast and window is None else window
        data = self.dists = self.__dtwDists(window)
        cls.fit(data)
        self.cluster_centers = [self.data[i] for i in cls.cluster_centers_indices_]
        self.clusters = collections.defaultdict(list)
        self.ass = list(cls.labels_)
        for ind, label in enumerate(cls.labels_):
            self.clusters[label].append(ind)
        self.n_clusters = len(self.clusters)

    def __saxDists(self):
        N = len(self.indexes)
        dist = np.zeros((N,N))
        for i in xrange(N):
            for j in xrange(i+1, N):
                x, y = self.indexes.keys()[i], self.indexes.keys()[j]
                dist[i, j] = dist[j, i] = st.dist(self.lookup, x, y)
        return dist

    def __dtwDists(self, window=None):
        N = self.N
        dist = np.zeros((N, N))
        for i in xrange(N):
            for j in xrange(i + 1, N):
                x, y = self.data[i, :], self.data[j, :]
                dist[i, j] = dist[j, i] = dtwf(x, y, euclidean)[0] if window else dtws(x, y, euclidean)[0]
        return dist

    def getDist(self,n1, n2):
        x, y = self.cluster_sax[n1], self.cluster_sax[n2]
        return self.dists[frozenset({x,y})]

