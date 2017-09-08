import numpy as np
import sdata as sd


def simulate(func, w, n, L):
    f = np.vectorize(func)
    xx = np.asarray([[x for x in range(-L,L +1)] for _ in xrange(n)])
    d = xx + np.random.random((n, 2*L+1)) * 2 * w - w
    y = sd.znorm(f(d)) + np.random.random((n, 2*L + 1))
    return y