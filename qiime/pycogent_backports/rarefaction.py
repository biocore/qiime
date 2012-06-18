#!/usr/bin/env python
from numpy import concatenate, repeat, array, zeros, histogram, arange, uint, zeros
from numpy.random import permutation, randint

"""Given array of objects (counts or indices), perform rarefaction analyses."""

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

def subsample(counts, n):
    """Subsamples new vector from vector of orig items.
    
    Returns all items if requested sample is larger than number of items.
    """
    if counts.sum() <= n:
        return counts
    nz = counts.nonzero()[0]
    unpacked = concatenate([repeat(array(i,), counts[i]) for i in nz])
    permuted = permutation(unpacked)[:n]
    result = zeros(len(counts))
    for p in permuted:
        result[p] += 1
    return result

def subsample_freq_dist_nonzero(counts, n, dtype=uint):
    """Subsamples new vector from vector of orig items.

    Returns all items if requested sample is larger than number of items.

    This version uses the cumsum/frequency distribution method.
    """
    if counts.sum() <= n:
        return counts
    result = zeros(len(counts), dtype=dtype)
    nz = counts.nonzero()[0]
    compressed = counts.take(nz)
    sums = compressed.cumsum()
    total = sums[-1]
    del compressed
    curr = n
    while curr:
        pick = randint(0, total)
        #print pick, sums, sums.searchsorted(pick), '\n'
        index = sums.searchsorted(pick,side='right')
        result[nz[index]] += 1
        sums[index:] -= 1
        curr -= 1
        total -= 1
    return result

def naive_histogram(vals, max_val=None, result=None):
    """Naive histogram for performance testing vs. numpy's.
    
    Apparently numpy's is 3x faster (non-cumulative) for larger step sizes
    (e.g. 1000) and 10x slower for small step sizes (e.g. 1), so will use
    logic to switch over depending on conditions.
    """
    if max_val is None:
       max_val = vals.max()
    if result is None:
        result = zeros(max_val+1, dtype=int)
    for v in vals:
        result[v] += 1
    return result

def wrap_numpy_histogram(max_val):
    """return convenience wrapper for numpy histogram"""
    bins = arange(max_val+2, dtype = int)   #+1 for length, +1 for leading 0
    def f(vals, max_val='ignored'): return histogram(vals, bins)[0]
    return f

def rarefaction(data, start=0, stop=None, stride=1, histogram_f=None, \
    permutation_f=permutation, is_counts=True):
    """Yields successive subsamples as vectors from vector of orig items.
    
    data can either be array of counts or array of observations. Default is
    to assume counts; set is_counts to False if this is not the case for your
    input.

    Returns all items if requested sample is larger than number of items.

    WARNING: each successive result is written into the same object (for
    convenience) so if you want the actual vectors for each rarefaction you
    need to do something like res = [r.copy() for r in rarefaction(params)].
    """
    if is_counts:   #need to transform data into indices
        nz = array(data).nonzero()[0]
        indices = concatenate([repeat(array(i,), data[i]) for i in nz])
    else:
        indices = array(data)

    if stop is None:
        stop = len(indices)
    if not stride:
        stride = 1  #avoid zero or None as stride
    max_val=indices.max() 
    if histogram_f is None:
        if stride < 100:
            histogram_f = naive_histogram
        else:
            histogram_f = wrap_numpy_histogram(max_val)
    permuted = permutation_f(indices)
    result = zeros(max_val+1, dtype=int)
    while start < stop:
        curr_slice = permuted[start:start+stride]
        result += histogram_f(curr_slice, max_val=max_val)
        yield result
        start += stride

