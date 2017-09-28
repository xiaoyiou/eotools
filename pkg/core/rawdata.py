"""
This module contains the classes for raw data containers
"""
from pybedtools import BedTool
import gc
import pandas as pd
from concurrent import futures
import cPickle as pickle
import numpy as np
from scipy.stats import poisson as pois


class Experiment(object):
    """
    This stores the processed data using genome and
    windows data
    """
    def __init__(self, paths=[], names=[]):
        assert len(paths) == len(names)
        self.paths = paths
        self.names = names
        self.raw = None
        self.binary = None
        self.genes = None

    def addTrack(self, path, name):
        self.paths.append(path)
        self.names.append(name)

    def preprocess(self, chrlenPath, genomePath, w=100, upStream=1000,
            downStream=1000, overlap=0.5, method='mean', col=4, type='bed', n_workers=4):


        assert upStream % w == 0 and downStream % w == 0
        window = BedTool().window_maker(g=chrlenPath, w=w)
        genes = BedTool(genomePath).to_dataframe()
        genes = genes[genes['feature'] == 'gene'][['seqname', 'start', 'strand', 'attributes']]
        genes['attributes'] = genes['attributes'].apply(lambda x: x[x.find('=')+1:x.find(';')])
        genes['start'] = genes['start'].apply(lambda x: x-upStream)
        genes['end'] = genes['start'] + upStream + downStream
        genes = genes[['seqname', 'start', 'end', 'attributes', 'strand']]
        genes.columns = ['chrom', 'start', 'end', 'ID', 'strand']
        genes = genes[genes.start >= 0]
        genes.chrom = genes.chrom.apply(lambda x:x[0].lower() + x[1:])
        atlas = BedTool.from_dataframe(genes[['chrom', 'start', 'end', 'ID']]).sort()
        genes = genes.set_index(['ID'])

        def worker(atlas, window, path, genes, col, method, overlap, type):
            p = BedTool(path).sort()
            a = None
            if type == 'bed':
                a = window.map(p, c=1, o='count', F=overlap)
            elif type == 'sigbed':
                a = window.map(p, o=method, c=col, F=overlap)
            tmp = atlas.intersect(a, loj=True, wa=True, wb=True).to_dataframe()
            grps = tmp.groupby(['name'])
            data = []
            for ind in genes.index:
                row = grps.get_group(ind)['thickEnd'].tolist()
                data.append(row if genes.ix[ind].strand == '+' else row[::-1])
            return pd.DataFrame(data, index=genes.index.tolist())
        self.raw = {}
        with futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
            jobs = {}
            for i, path in enumerate(self.paths):
                job = executor.submit(worker, atlas, window, path, genes, col, method, overlap, type)
                jobs[job] = self.names[i]

            for job in futures.as_completed(jobs):
                self.raw[jobs[job]] = job.result().dropna()
                if self.genes == None:
                    self.genes = self.raw[jobs[job]].index.tolist()


    def export(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)
            f.close()
        gc.collect()


    def get_binary(self, power=4):
        if self.raw is None or len(self.raw) == 0:
            return
        self.binary = pd.DataFrame(np.zeros((len(self.genes), len(self.names)),dtype=np.int),
                                   index = self.genes,
                                   columns = self.names
                                   )
        print self.raw
        for mod in self.raw:
            print "debugging!!!", mod
            print self.raw[mod].values.mean()
            x = pois(self.raw[mod].values.mean())
            thresh = x.isf(10**(-power))
            self.binary.ix[(self.raw[mod] >= thresh).any(axis=1), [mod]] = 1
        self.binary.columns = range(len(self.names))
