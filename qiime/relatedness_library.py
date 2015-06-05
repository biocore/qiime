#!/usr/bin/env python
# File created on 02 May 2012
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Van Treuren"
__email__ = "wdwvt1@gmail.com"

from numpy.random import shuffle
from numpy import std, mean, arange, eye
from numpy.ma import masked_array


"""The calculations for MPD (mean phylogenetic distance), MNTD (mean nearest
taxon distance), NRI (net relatedness index), and NTI (nearest taxon index) are
based on the formulas available in the Phylocom 4.2 manual and the Webb et al.
2002 paper in Phylogenies and Comunnity Ecology. The manual is
currently available at http://phylodiversity.net/phylocom/phylocom_manual.pdf,
and the paper is available at
http://www.annualreviews.org/doi/pdf/10.1146/annurev.ecolsys.33.010802.150448.

The null model that is used by the random_mpd and random_mntd methods is null
model 2 (pg 16). Null model 2 specifies that for random draws from the total
phylogeny, any taxa may be taken regardless of whether or not it occurred in
one of the samples.


Comparisons against phylocom 4.2 shows significant agreement. If you would
like to test, make the sample file have the folowing form:
clump1    sp1
clump1    sp2
clump1    sp4
clump1    sp7

and the phylo file have the following tree string:
(((sp1:.06,sp2:.1)A:.031,(sp3:.001,sp4:.01)B:.2)AB:.4,((sp5:.03,sp6:.02)C:.13,(sp7:.01,sp8:.005)D:.1)CD:.3)root;

tr = DndParser(ts)
distmat = tr.tipToTipDistances()
distmat =
(array([[ 0.   ,  0.16 ,  0.292,  0.301,  0.951,  0.941,  0.901,  0.896],
       [ 0.16 ,  0.   ,  0.332,  0.341,  0.991,  0.981,  0.941,  0.936],
       [ 0.292,  0.332,  0.   ,  0.011,  1.061,  1.051,  1.011,  1.006],
       [ 0.301,  0.341,  0.011,  0.   ,  1.07 ,  1.06 ,  1.02 ,  1.015],
       [ 0.951,  0.991,  1.061,  1.07 ,  0.   ,  0.05 ,  0.27 ,  0.265],
       [ 0.941,  0.981,  1.051,  1.06 ,  0.05 ,  0.   ,  0.26 ,  0.255],
       [ 0.901,  0.941,  1.011,  1.02 ,  0.27 ,  0.26 ,  0.   ,  0.015],
       [ 0.896,  0.936,  1.006,  1.015,  0.265,  0.255,  0.015,  0.   ]]),
100 tests (87 for phylocom) using 1000 iters each:
nri
phylocom - mean = .4481, std = .05529
rel - mean = .4516, std = .048414

nti
phylocom - mean = -1.2081908, std = 0.03177133
rel - mean = -1.23526, std = 0.02437
"""

# used by both NRI and NTI


def reduce_mtx(distmat, indices):
    """Returns rows,cols of distmat where rows,cols=indices."""
    return distmat.take(indices, 0).take(indices, 1)

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
    if abs(sd_x_n) < .00001:
        raise ValueError('The standard deviation of the means of the random' +
                         ' draws from the distance matrix was less than .00001. This is' +
                         ' likely do to a phylogeny with distances to all tips equal.' +
                         ' This phylogeny is not suitable for NRI/NTI analysis.')
    return -1. * ((mn_x_obs - mn_x_n) / sd_x_n)


def mpd(distmat):
    """Return mean of pairwise dists, assumes distmat is symmetric,hollow."""
    return distmat.sum() / (distmat.size - distmat.shape[0])


def random_mpd(distmat, n, iters):
    """Calc mean,std of mean of iters # of rand nxn distmats drawn from distmat.
    Notes:
     Calculate the std of the means of the distances in the nxn mat excluding
     the main diagonal and below since d(A,A)=0 and distmat is symmetric.
     The forumula from Webb 2002 seems to calculate the standard deviation of
     the distances but based on tests of Phylocom and the Phylocom manual,
     Phylocom computes the standard deviations of the means."""
    means = []
    indices = arange(distmat.shape[0])  # square so rows=cols
    for i in range(iters):
        shuffle(indices)  # shuffling indices after its been shuffled is not
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
    if abs(sd_y_n) < .00001:
        raise ValueError('The standard deviation of the means of the random' +
                         ' draws from the distance matrix was less than .00001. This is' +
                         ' likely do to a phylogeny with distances to all tips equal.' +
                         ' This phylogeny is not suitable for NRI/NTI analysis.')
    return -1. * ((mn_y_obs - mn_y_n) / sd_y_n)


def mntd(distmat):
    """Find mean of row mins in hollow, symmetric distmat excluding main diag."""
    return masked_array(distmat, eye(distmat.shape[0])).min(0).mean()


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
    indices = arange(distmat.shape[0])  # square so rows=cols
    for i in range(iters):
        shuffle(indices)  # shuffling indices after its been shuffled is not
        # mathematically different than shuffling fresh arange(n)
        means.append(mntd(reduce_mtx(distmat, indices[:n])))
    return mean(means), std(means)
