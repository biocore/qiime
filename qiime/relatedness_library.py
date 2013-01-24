#!/usr/bin/env python
# File created on 02 May 2012
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "William Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"
 
from numpy.random import shuffle, random
from numpy import std, mean, array, allclose, arange, eye
from copy import deepcopy
import numpy.ma as ma

from numpy import triu_indices

"""The calculations for MPD (mean phylogenetic distance), MNTD (mean nearest
taxon distance), NRI (net relatedness index), and NTI (nearest taxon index) are
based on the formulas available in the Phylocom 4.2 manual. That manual is 
currently available at http://phylodiversity.net/phylocom/phylocom_manual.pdf.

The null model that is used by the random_mpd and random_mntd methods is null
model 2 (pg 16). Null model 2 specifies that for random draws from the total
phylogeny, any taxa may be taken regardless of whether or not it occurred in 
one of the samples.
"""



# used by both NRI and NTI

def reduce_mtx(distmat, indices):
    """Returns rows,cols of distmat where rows,cols=indices."""
    return distmat.take(indices,0).take(indices,1)

# NRI

def nri(distmat, marginals, group, iters):
    """Calculate the NRI of the selected group.
    Notes:
     distmat - distance matrix between taxa.
     marginals - list of ids of the marginals of the distmat, i.e. the taxa 
     names. only need to provide row xor col marginals since its a dist mat.
     group - list of ids of the group that you want to calculate the nri of. 
     iters - number of iterations to use. 1000 is suggested.
     MPD (mean phylogenetic distance) is calculated by mpd. 
    """
    group_marginals = [marginals.index(i) for i in group] 
    mn_x_obs = mpd(reduce_mtx(distmat, group_marginals))
    mn_x_n, sd_x_n = random_mpd(distmat, len(group_marginals), iters)
    if abs(0.0 - sd_x_n)<.00001:
        raise ValueError('The standard deviation of the means of the random'+\
            ' draws from the distance matrix was less than .00001. This is'+\
            ' likely do to a phylogeny with distances to all tips equal.'+\
            ' This phylogeny is not suitable for NRI/NTI analysis.')
    return -1.*((mn_x_obs-mn_x_n)/sd_x_n)

def mpd(distmat):
    """Return mean distmat, assumes distmat is symmetric and hollow."""
    return distmat.sum()/(distmat.size-distmat.shape[0])

def random_mpd(distmat, n, iters):
    """Calc mean,std of mean of iters # of rand nxn distmats drawn from distmat.
    Notes:
     Calculate the std of the means of the distances in the nxn mat excluding 
     the main diagonal and below since d(A,A)=0 and distmat is symmetric.
     The forumula from Webb 2002 seems to calculate the standard deviation of 
     the distances but based on tests of Phylocom and the Phylocom manual, 
     Phylocom computes the standard deviations of the means."""
    means = []
    indices = arange(distmat.shape[0]) #square so rows=cols
    for i in range(iters):
        shuffle(indices) #shuffling indices after its been shuffled is not 
        # mathematically different than shuffling fresh arange(n) 
        means.append(mpd(reduce_mtx(distmat, indices[:n])))
    return mean(means), std(means)

# NTI

def nti(distmat, marginals, group, iters):
    """Calculates the NTI of the selected group.
    Notes:
     distmat - distance matrix between taxa.
     marginals - list of ids of the marginals of the distmat, i.e. the taxa 
     names. only need to provide row xor col marginals since its a dist mat.
     group - list of ids of the group that you want to calculate the nti of. 
     iters - number of iterations to use. 1000 is suggested.
     MNTD (mean nearest taxon distance) is calculated by mntd. 
    """
    group_marginals = [marginals.index(i) for i in group] 
    mn_y_obs = mntd(reduce_mtx(distmat, group_marginals))
    mn_y_n, sd_y_n = random_mntd(distmat, len(group_marginals), iters)
    if abs(0.0 - sd_y_n)<.00001:
        raise ValueError('The standard deviation of the means of the random'+\
            ' draws from the distance matrix was less than .00001. This is'+\
            ' likely do to a phylogeny with distances to all tips equal.'+\
            ' This phylogeny is not suitable for NRI/NTI analysis.')
    return -1.*((mn_y_obs-mn_y_n)/sd_y_n)

def mntd(distmat):
    """Find mean of row mins in hollow, symmetric distmat excluding main diag."""
    return ma.masked_array(distmat, eye(distmat.shape[0])).min(0).mean()

def random_mntd(distmat, n, iters):
    """Calc mean,std of mntd of iters # of rand nxn mtx's drawn from distmat.
    Notes:
     Calculate the std of the means of minimums of the distances in the nxn mat
     excluding the main diagonal since d(A,A)=0 (mtx is hollow).
     The forumula from Webb 2002 seems to calculate the standard deviation of 
     the distances but based on tests of Phylocom and the Phylocom manual, 
     Phylocom computes the standard deviations of the means.
     """
    means = []
    indices = arange(distmat.shape[0]) #square so rows=cols
    for i in range(iters):
        shuffle(indices) #shuffling indices after its been shuffled is not 
        # mathematically different than shuffling fresh arange(n) 
        means.append(mntd(reduce_mtx(distmat, indices[:n])))
    return mean(means), std(means)


