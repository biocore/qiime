#!/usr/bin/env python
from numpy import concatenate, repeat, array, zeros, histogram, arange, uint, zeros
from numpy.random import permutation, randint, sample
from random import Random, _ceil, _log

"""Given array of objects (counts or indices), perform rarefaction analyses."""

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"


class MyRandom(Random):
    """Adding a method to sample from array only"""
    def sample_array(self, population, k):
        """Chooses k unique random elements from a population sequence.

        Returns a new list containing elements from the population while
        leaving the original population unchanged.  The resulting list is
        in selection order so that all sub-slices will also be valid random
        samples.  This allows raffle winners (the sample) to be partitioned
        into grand prize and second place winners (the subslices).

        Members of the population need not be hashable or unique.  If the
        population contains repeats, then each occurrence is a possible
        selection in the sample.

        To choose a sample in a range of integers, use xrange as an argument.
        This is especially fast and space efficient for sampling from a
        large population:   sample(xrange(10000000), 60)
        """

        # Sampling without replacement entails tracking either potential
        # selections (the pool) in a list or previous selections in a set.

        # When the number of selections is small compared to the
        # population, then tracking selections is efficient, requiring
        # only a small set and an occasional reselection.  For
        # a larger number of selections, the pool tracking method is
        # preferred since the list takes less space than the
        # set and it doesn't suffer from frequent reselections.

        n = len(population)
        if not 0 <= k <= n:
            raise ValueError("sample larger than population")
        random = self.random
        _int = int
        result = zeros(k)
        setsize = 21        # size of a small set minus size of an empty list
        if k > 5:
            setsize += 4 ** _ceil(_log(k * 3, 4)) # table size for big sets
        if n <= setsize or hasattr(population, "keys"):
            # An n-length list is smaller than a k-length set, or this is a
            # mapping type so the other algorithm wouldn't work.
            pool = array(list(population))
            for i in xrange(k):         # invariant:  non-selected at [0,n-i)
                j = _int(random() * (n-i))
                result[i] = pool[j]
                pool[j] = pool[n-i-1]   # move non-selected item into vacancy
        else:
            try:
                selected = set()
                selected_add = selected.add
                for i in xrange(k):
                    j = _int(random() * n)
                    while j in selected:
                        j = _int(random() * n)
                    selected_add(j)
                    result[i] = population[j]
            except (TypeError, KeyError):   # handle (at least) sets
                if isinstance(population, list):
                    raise
                return self.sample_array(tuple(population), k)
        return result

_inst = MyRandom()
sample = _inst.sample_array


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

def subsample_random(counts, n, dtype=uint):
    """Subsamples new vector from vector of orig items.

    Returns all items if requested sample is larger than number of items.

    This version uses random.sample.
    """
    if counts.sum() <= n:
        return counts
    nz = counts.nonzero()[0]
    unpacked = concatenate([repeat(array(i,), counts[i]) for i in nz])
    permuted = sample(unpacked, n)
    result = zeros(len(counts),dtype=dtype)
    for p in permuted:
        result[p] += 1
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

