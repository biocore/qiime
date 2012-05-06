#!/usr/bin/env python
from __future__ import division

__author__ = "Rob Knight, Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"


"""contains metrics for use with beta_diversity.py
most metrics are from cogent.math.distance_transform.py,
but some need wrappers to look like f(data, taxon_names, tree)-> dist_mtx
"""
import cogent.maths.unifrac.fast_tree as fast_tree
 # (unifrac, unnormalized_unifrac,
 #    G, unnormalized_G, weighted_unifrac)
from cogent.maths.unifrac.fast_unifrac import fast_unifrac, fast_unifrac_one_sample
from qiime.parse import make_envs_dict
import numpy
import warnings

# add pycogent's bray-curtis as bray-curtis-faith, note alternate formula exists
from cogent.maths.distance_transform import dist_bray_curtis as dist_bray_curtis_faith

def dist_bray_curtis_magurran(datamtx, strict=True):
    """ returns bray curtis distance (quantitative sorensen) btw rows
    
    dist(a,b) = 2*sum on i( min( a_i, b_i)) / sum on i( (a_i + b_i) )
    
    see for example:
    Magurran 2004
    Bray 1957

    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * if strict==True, raises ValueError if any of the input data is negative,
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a matrix with nonnegative 
    entries.  If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not numpy.all(numpy.isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if numpy.any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if numpy.rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = numpy.shape(datamtx)
    else:
        try:
            numrows, numcols = numpy.shape(datamtx)
        except ValueError:
            return numpy.zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return numpy.zeros((0,0),'d')

    dists = numpy.zeros((numrows,numrows),'d')
    for i in range(numrows):
        r1 = datamtx[i,:]
        r1sum = r1.sum()
        for j in range(i):
            r2 = datamtx[j,:]
            r2sum = r2.sum()
            minvals = numpy.min([r1,r2],axis=0)

            if (r1sum + r2sum) == 0:
                dists[i][j] = dists[j][i] = 0.0
            else:
                dissim = 1 - ( (2*minvals.sum()) / (r1sum + r2sum) )
                dists[i][j] = dists[j][i] = dissim
    return dists


def make_unifrac_metric(weighted, metric, is_symmetric):
    """Make a unifrac-like metric.

    Parameters:
    sample_names: list of unique strings
    weighted: bool, if True does the necessary root to tip calcs.
    metric: f(branch_lengths, i, j) -> distance
    is_symmetric: saves calc time if metric is symmetric.
    kwargs passed to fast_unifrac
    """
    def result(data, taxon_names, tree, sample_names, **kwargs):
        """ wraps the fast_unifrac fn to return just a matrix, in correct order
        
            sample_names: list of unique strings
        """

        envs = make_envs_dict(data, sample_names, taxon_names)
        unifrac_res = fast_unifrac(tree, envs, weighted=weighted, metric=metric,
            is_symmetric=is_symmetric, modes=["distance_matrix"],**kwargs)
        dist_mtx = _reorder_unifrac_res(unifrac_res['distance_matrix'],
            sample_names)
        return dist_mtx
    return result

# these should start with dist_ to be discoverable by beta_diversity.py
# unweighted full tree => keep the full tree relating all samples.  
# Compute how much branch
# length is present in one sample but not (both samples OR NEITHER
# SAMPLE).  Divide by total branch length of full tree.  
# G is asymmetric unifrac
dist_unweighted_unifrac = make_unifrac_metric(False, fast_tree.unifrac, True)
dist_unifrac = dist_unweighted_unifrac # default unifrac is just unifrac
dist_unweighted_unifrac_full_tree = make_unifrac_metric(False, 
    fast_tree.unnormalized_unifrac, True)
dist_weighted_unifrac = make_unifrac_metric(True, 
    fast_tree.weighted_unifrac, True)
dist_weighted_normalized_unifrac = make_unifrac_metric('correct',
    fast_tree.weighted_unifrac, True)
dist_unifrac_g = make_unifrac_metric(False, fast_tree.G, False)
dist_unifrac_g_full_tree = make_unifrac_metric(False, 
    fast_tree.unnormalized_G, False)

def make_unifrac_row_metric(weighted, metric, is_symmetric):
    """Make a unifrac-like metric, for only one row of the dissm mtx

    Parameters:
    sample_names: list of unique strings
    weighted: bool, if True does the necessary root to tip calcs.
    metric: f(branch_lengths, i, j) -> distance
    is_symmetric: ignored
    sample_name: of the sample corresponding to the row of the dissim mtx
    """
    def result(data, taxon_names, tree, sample_names, one_sample_name,**kwargs):
        """ wraps the fast_unifrac fn to return just a matrix, in correct order

            sample_names: list of unique strings
        """
        envs = make_envs_dict(data, sample_names, taxon_names)
        unifrac_res = fast_unifrac_one_sample(one_sample_name,
            tree, envs, weighted=weighted, metric=metric,**kwargs)
        dist_mtx = _reorder_unifrac_res_one_sample(unifrac_res,
            sample_names)
        return dist_mtx
    return result

one_sample_unweighted_unifrac = make_unifrac_row_metric(False, fast_tree.unifrac, True)
one_sample_unifrac = one_sample_unweighted_unifrac # default unifrac is just unifrac
one_sample_unweighted_unifrac_full_tree = make_unifrac_row_metric(False, 
    fast_tree.unnormalized_unifrac, True)
one_sample_weighted_unifrac = make_unifrac_row_metric(True, 
    fast_tree.weighted_unifrac, True)
one_sample_weighted_normalized_unifrac = make_unifrac_row_metric('correct',
    fast_tree.weighted_unifrac, True)
one_sample_unifrac_g = make_unifrac_row_metric(False, fast_tree.G, False)
one_sample_unifrac_g_full_tree = make_unifrac_row_metric(False, 
    fast_tree.unnormalized_G, False)

def _reorder_unifrac_res(unifrac_res, sample_names_in_desired_order):
    """ reorder unifrac result
    
    unifrac res is distmtx,sample_names.  sample names not in unifrac's
    sample names (not in tree, all zeros in otu table(?)) will be included, 
    with a user warning.
    
    """
    sample_names = sample_names_in_desired_order
    unifrac_dist_mtx = unifrac_res[0]
    unifrac_sample_names = unifrac_res[1]
    if unifrac_sample_names == sample_names:
        dist_mtx = unifrac_dist_mtx
    else:
        dist_mtx = numpy.zeros((len(sample_names), len(sample_names)))
        
        for i, sam_i in enumerate(sample_names):
            
            # make dist zero if both absent, else dist=1. dist to self is 0
            if sam_i not in unifrac_sample_names:
                warnings.warn('unifrac had no information for sample ' +\
                 sam_i + ". Distances involving that sample aren't meaningful")
                for j, sam_j in enumerate(sample_names):
                    if sam_j not in unifrac_sample_names:
                        dist_mtx[i,j] = 0.0
                    else:
                        dist_mtx[i,j] = 1.0
                        
            # sam_i is present, so get unifrac dist
            else:
                unifrac_i =  unifrac_sample_names.index(sam_i)
                
                for j, sam_j in enumerate(sample_names):
                    if sam_j not in unifrac_sample_names:
                        dist_mtx[i,j] = 1.0
                    else:
                        unifrac_j =  unifrac_sample_names.index(sam_j)
                        dist_mtx[i,j] = unifrac_dist_mtx[unifrac_i, unifrac_j]
    return dist_mtx


def _reorder_unifrac_res_one_sample(unifrac_res, sample_names_in_desired_order):
    """ reorder unifrac result
    
    unifrac res is distmtx,sample_names.  sample names not in unifrac's
    sample names (not in tree, all zeros in otu table(?)) will be included, 
    with a user warning.
    
    """
    sample_names = sample_names_in_desired_order
    unifrac_dist_arry = unifrac_res[0]
    unifrac_sample_names = unifrac_res[1]
    if unifrac_sample_names == sample_names:
        dist_arry = unifrac_dist_arry
    else:
        dist_arry= numpy.zeros(len(sample_names))
        
        for i, sam_i in enumerate(sample_names):
            
            # make dist zero if both absent, else dist=1. dist to self is 0
            if sam_i not in unifrac_sample_names:
                warnings.warn('unifrac had no information for sample ' +\
                 sam_i + ". Distances involving that sample aren't meaningful")
                dist_arry[i] = 1.0
                        
            # sam_i is present, so get unifrac dist
            else:
                unifrac_i =  unifrac_sample_names.index(sam_i)
                dist_arry[i] = unifrac_dist_arry[unifrac_i]
    return dist_arry
