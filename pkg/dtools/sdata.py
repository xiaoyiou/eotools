"""
This module contains the necessary tools for preprocesing series like datasets. (s-data)
In all the array-like data, the rows are actual entries of data points where each row represent a gene
. Approximation (shrinking the length of each vector) methods using different; the returned value is still
a numpy array
. Z-transformation to normalize the data;
. SAX symbolic transformation;
"""
import numpy as np
import pandas as pd
import scipy.stats.mstats as stats
from scipy.stats import norm
import pywt
from scipy.spatial.distance import euclidean
# Private methods:


def __paaProx(ts, n_pieces):
    """
    ts: the columns of which are time series represented by e.g. np.array
    n_pieces: M equally sized piecies into which the original ts is splitted
    """
    ts = np.asarray(ts, dtype=np.float)
    splitted = np.array_split(ts, n_pieces, axis=1) # along columns as we want
    return np.asarray(map(lambda xs: xs.mean(axis=1), splitted)).T


def __pmaProx(ts, n_pieces):
    """

    :param ts: A row-vector based series data,
    :param n_pieces: window size
    :return: the shortened data where each value is the maximum of window
    """
    splitted = np.array_split(ts, n_pieces, axis=1)  # along columns as we want
    return np.asarray(map(lambda xs: xs.max(axis=1), splitted)).T


def __dwtProx(ts, level, wavelet):
    w = pywt.Wavelet(wavelet)
    return pywt.waverec(pywt.wavedec(ts, w, mode='constant')[:level], w)


def znorm (ts):
    return stats.zscore(ts, axis=1)


def approximate(ts, n, method='dwt', wavelet='db1'):
    """
    :param ts: Input row-vector based data
    :param n: number of points remaining after approximation. Note that if method = dwt, then the returned
    new dimension will be the largest 2^k <= n. The largests power of 2 less than or equal to n.
    :param method: 'paa, pma, dwt'
    :return: the processed data
    """
    res = None
    if method == 'dwt':
        level = n.bit_length() - 1
        res = __dwtProx(ts, level, wavelet)
    elif method == 'paa':
        res = __paaProx(ts, n)
    elif method == 'pma':
        res = __pmaProx(ts, n)
    return res


def dist(dist, s1, s2 , p=2):
    assert len(s1) == len(s2)
    if dist is None:
        return euclidean(s1, s2)
    total = 0
    for i in xrange(len(s1)):
        total += dist.loc[s1[i], s2[i]] ** p
    return total ** (1.0/p)


def sax_transform(ts, alphabet, thrholds=None):
    """
    Before using this, the ts should be compressed (suing wavelet or paa)
    and then z transformed to ensure equal frequency assignments
    Return: (The sax representations matrix, distance matrix)
    """
    alphabet_sz = len(alphabet)
    if thrholds is None:
        thrholds = norm.ppf(np.linspace(1./alphabet_sz,
                                    1-1./alphabet_sz,
                                   alphabet_sz-1))

    def translate(ts_values):
        return np.asarray([(alphabet[0] if ts_value < thrholds[0]
                else (alphabet[-1] if ts_value > thrholds[-1]
                      else alphabet[np.where(thrholds <= ts_value)[0][-1]+1]))
                           for ts_value in ts_values])
    dist = np.zeros((alphabet_sz, alphabet_sz))

    for i in xrange(alphabet_sz - 1):
        for j in xrange(i + 1, alphabet_sz):
            if abs(i - j) <= 1:
                continue
            dist[i,j] = thrholds[j-1] - thrholds[i]
            dist[j,i] = dist[i,j]

    return (np.apply_along_axis(translate, 0, ts),\
            pd.DataFrame(dist, index=list(alphabet), columns=list(alphabet)),\
            thrholds)
