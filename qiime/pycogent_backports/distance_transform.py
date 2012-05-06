#!/usr/bin/env python
""" matrix based distance metrics, and related coordinate transforms
     
functions to compute distance matrices row by row from abundance matrices,
typically samples (rows) vs. species/OTU's (cols)

DISTANCE FUNCTIONS
For distance functions, the API resembles the following (but see function 
docstring for specifics):
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros *typically* returns 0 distance between them
    * negative values are only allowed for some distance metrics,
    in these cases if strict==True, negative input values return a ValueError, 
    and if strict==False, errors or misleading return values may result
    * functions prefaced with "binary" consider only presense/absense in
    input data (qualitative rather than quantitative)

TRANSFORM FUNCTIONS
* For transform functions, very little error checking exists.  0/0 evals
in transform fomulas will throw errors, and negative data will return 
spurious results or throw errors
* The transform functions are as described in
    Legendre, P. and E. Gallagher. 2001.  Ecologically meaningful
    transformations for ordination of species data.  Oecologia: 129: 271-280.
These and allow the use
of ordination methods such as PCA and RDA, which are Euclidean-based, 
for the analysis of community data, while circumventing the problems associated
with the Euclidean distance. The matrix that is returned still has samples as
rows and species as columns, but the values are transformed so that when
programs such as PCA calculate euclidean distances on the matrix, chord,
chisquare, 'species profile', or hellinger distances will result.

EXAMPLE USAGE: 
    >from distance_transform import dist_euclidean
    >from numpy import array

    >abundance_data = array([[1, 3],
                            [5, 2],
                            [0.1, 22]],'d')

    >dists = dist_euclidean(abundance_data)

    >print dists
    
    array([[  0.        ,   4.12310563,  19.02130385],
           [  4.12310563,   0.        ,  20.5915031 ],
           [ 19.02130385,  20.5915031 ,   0.        ]])

    
"""
from __future__ import division
import numpy
from numpy import (array, zeros, logical_and, logical_or, logical_xor, where,
    mean, std, argsort, take, ravel, logical_not, shape, sqrt, abs, 
    sum, square, asmatrix, asarray, multiply, min, rank, any, all, isfinite,
    nonzero, nan_to_num, geterr, seterr, isnan)
# any, all from numpy override built in any, all, preventing:
# ValueError: The truth value of an array with more than one element is 
# ambiguous. Use a.any() or a.all()

from numpy.linalg import norm

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Rob Knight", "Micah Hamady", "Justin Kuczynski",
                    "Zongzhi Liu", "Catherine Lozupone", 
                    "Antonio Gonzalez Pena", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

def _rankdata(a):
    """ Ranks the data in a, dealing with ties appropritely.  First ravels 
    a.  Adapted from Gary Perlman's |Stat ranksort.
    private helper function

    Returns: array of length equal to a, containing rank scores
    """
    a = ravel(a)
    n = len(a)
    ivec = argsort(a)
    svec = take(a, ivec)
    sumranks = dupcount = 0
    newarray = zeros(n,'d')
    for i in range(n):
        sumranks = sumranks + i
        dupcount = dupcount + 1
        if i==n-1 or svec[i] <> svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in range(i-dupcount+1,i+1):
                newarray[ivec[j]] = averank
                sumranks = dupcount = 0
    return newarray

def trans_chord(m):
    """perform a chord distance transformation on the rows of m

    transforms m to m' so that the euclidean dist between the rows of m' equals
    the chord dist between the rows of m.
    Ref:
    Legendre, P. and E. Gallagher. 2001.  Ecologically meaningful
    transformations for ordination of species data.  Oecologia: 129: 271-280.
    """
    m = asmatrix(m)
    row_norms = sqrt(sum(square(m), axis=1))
    result = m / row_norms
    return result

def trans_chisq(m):
    """perform a chi squared distance transformation on the rows of m

    transforms m to m' so that the euclidean dist between the rows of m' equals
    the chi squared dist between the rows of m.    
    Ref:
    Legendre, P. and E. Gallagher. 2001.  Ecologically meaningful
    transformations for ordination of species data.  Oecologia: 129: 271-280.
    """
    m = asmatrix(m)
    grand_sum, row_sums, col_sums = m.sum(), m.sum(1), m.sum(0)
    result = m * sqrt(grand_sum)
    result /= row_sums
    result /= sqrt(col_sums)
    return result

def trans_specprof(m):
    """perform a species profile distance transformation on the rows of m

    transforms m to m' so that the euclidean dist between the rows of m' equals
    the species profile dist between the rows of m.
    Ref:
    Legendre, P. and E. Gallagher. 2001.  Ecologically meaningful
    transformations for ordination of species data.  Oecologia: 129: 271-280.
    """
    m = asmatrix(m)
    row_sums = sum(m, axis=1)
    result = m / row_sums
    return result

def trans_hellinger(m):
    """perform a hellinger distance transformation on the rows of m

    transforms m to m' so that the euclidean dist between the rows of m' equals
    the hellinger dist between the rows of m.
    Ref:
    Legendre, P. and E. Gallagher. 2001.  Ecologically meaningful
    transformations for ordination of species data.  Oecologia: 129: 271-280.
    """
    m = asmatrix(m)
    row_sums = sum(m, axis=1)
    result = sqrt(m / row_sums)
    return result



def dist_bray_curtis(datamtx, strict=True):
    """ returns bray curtis distance (normalized manhattan distance) btw rows
    
    dist(a,b) = manhattan distance / sum on i( (a_i + b_i) )
    
    see for example:
    Faith et al. 1987
    Compositional dissimilarity as a robust measure of ecological distance
    Vegitatio

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
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')

    dists = zeros((numrows,numrows),'d')
    for i in range(numrows):
        r1 = datamtx[i,:]
        for j in range(i):
            r2 = datamtx[j,:]
            abs_v = float(sum(abs(r1 - r2)))
            v = sum(r1 + r2)
            cur_d = 0.0 
            if v > 0:
                cur_d = abs_v/v

            dists[i][j] = dists[j][i] = cur_d
    return dists

dist_bray_curtis_faith = dist_bray_curtis

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

def dist_canberra(datamtx, strict=True):
    """returns a row-row canberra dist matrix
    
    see for example:
    Faith et al. 1987
    Compositional dissimilarity as a robust measure of ecological distance
    Vegitatio

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
    * chisq dist normalizes by column sums - empty columns (all zeros) are
    ignored here
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    oldstate = seterr(invalid='ignore',divide='ignore')
    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    for i in range(numrows):
        r1 = datamtx[i]
        for j in range(i):
            r2 = datamtx[j]
            dist = 0.0
            net = abs( r1 - r2 ) / (r1 + r2)

            net = nan_to_num(net)
            num_nonzeros = nonzero(net)[0].size
            dists[i,j] = dists[j,i] = nan_to_num(net.sum()/num_nonzeros)
            
    seterr(**oldstate)
    return dists

def dist_chisq(datamtx, strict=True):
    """returns a row-row chisq dist matrix

    see for example:
    Faith et al. 1987
    Compositional dissimilarity as a robust measure of ecological distance
    Vegitatio
    
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * an all zero row compared with a not all zero row returns a distance of 1
    * if strict==True, raises ValueError if any of the input data is negative,
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a matrix with nonnegative 
    entries.  If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    * chisq dist normalizes by column sums - empty columns (all zeros) are
    ignored here
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    sqrt_grand_sum = sqrt(sum(datamtx))
    rowsums, colsums = sum(datamtx, axis=1), sum(datamtx, axis=0)
    if not colsums.all():
        for i in range(len(colsums)):
            if colsums[i] == 0.0:
                colsums[i] = 1.0
    for i in range(numrows):
        r1 = datamtx[i]
        r1sum = rowsums[i]
        for j in range(i):
            r2 = datamtx[j]
            r2sum = rowsums[j]
            if r1sum == 0.0 or r2sum == 0.0:
                if r1sum == 0.0 and r2sum == 0.0:
                    dist = 0.0
                else: dist = 1.0
            else:
                dist = sqrt_grand_sum *\
                    sqrt(sum( multiply((1./colsums) ,
                    square(r1/r1sum - r2/r2sum)) ))
            dists[i,j] = dists[j,i] = dist
    return dists

def dist_chord(datamtx, strict=True):
    """returns a row-row chord dist matrix

    attributed to Orloci (with accent).  see Legendre 2001,
    ecologically meaningful... 
    
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * an all zero row compared with a not all zero row returns a distance of 1
    * if strict==True, raises ValueError if any of the input data is 
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a 2d matrix.  
    If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    for i in range(numrows):
        r1 = datamtx[i] # cache here
        r1norm = norm(r1)
        for j in range(i):
            r2 = datamtx[j]
            r2norm = norm(r2)
            if r1norm == 0.0 or r2norm == 0.0:
                if r1norm == 0.0 and r2norm == 0.0:
                    dist = 0.0
                else: dist = 1.0
            else:
                dist = norm(r1/r1norm - r2/r2norm)
            dists[i,j] = dists[j,i] = dist

    return dists

def dist_euclidean(datamtx, strict=True):
    """returns a row by row euclidean dist matrix
    
    returns the euclidean norm of row1 - row2 for all rows in datamtx
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * if strict==True, raises ValueError if any of the input data is 
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a 2d matrix.  
    If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    datamtx = asarray(datamtx, 'd')
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    for r in range(numrows):
        for c in range(r):
            dist = norm(datamtx[r] - datamtx[c])
            if isnan(dist):
                raise RuntimeError('ERROR: overflow when computing euclidean distance')
            dists[r,c] = dists[c,r] = dist

    return dists

def dist_gower(datamtx, strict=True):
    """returns a row-row gower dist matrix
    
    see for example, Faith et al., 1987
    
    
    * note that the comparison between any two rows is dependent on the entire
    data matrix, d_ij is a fn of all of datamtx, not just i,j
    * comparisons are between rows (samples)
    * any column containing identical data for all rows is ignored (this
    prevents a 0/0 error in the formula for gower distance
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * if strict==True, raises ValueError if any of the input data is
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a 2d matrix.  
    If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    coldiffs = datamtx.max(axis=0) - datamtx.min(axis=0)
    for i in range(numcols):
        if coldiffs[i] == 0.0:
            coldiffs[i] = 1.0 # numerator will be zero anyway

    for i in range(numrows):
        r1 = datamtx[i]
        for j in range(i):
            r2 = datamtx[j]
            rowdiff = r2 - r1
            dist = sum(abs(r1 - r2) / coldiffs)
            dists[i,j] = dists[j,i] = dist

    return dists

def dist_hellinger(datamtx, strict=True):
    """returns a row-row hellinger dist matrix
    
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * an all zero row compared with a not all zero row returns a distance of 1
    * if strict==True, raises ValueError if any of the input data is negative,
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a matrix with nonnegative 
    entries.  If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    for i in range(numrows):
        r1 = datamtx[i]
        r1sum = sum(r1)
        for j in range(i):
            r2 = datamtx[j]
            r2sum = sum(r2)
            if r1sum == 0.0 or r2sum == 0.0:
                if r1sum == 0.0 and r2sum == 0.0:
                    dist = 0.0
                else: dist = 1.0
            else:
                dist = norm(sqrt(r1/r1sum) - sqrt(r2/r2sum))
            dists[i,j] = dists[j,i] = dist

    return dists

def dist_kulczynski(datamtx, strict=True):
    """ calculates the kulczynski distances between rows of a matrix
    
    see for example Faith et al., composiitonal dissimilarity, 1987
    returns a distance of 1 between a row of zeros and a row with at least one
    nonzero element
    
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * an all zero row compared with a not all zero row returns a distance of 1
    * if strict==True, raises ValueError if any of the input data is negative,
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a matrix with nonnegative 
    entries.  If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    rowsums = datamtx.sum(axis=1)
    # rowsum: the sum of elements in a row
    # cache to avoid recalculating for each pair
    for i in range(numrows):
        irowsum = rowsums[i]
        r1 = datamtx[i]
        for j in range(i):
            r2 = datamtx[j]
            jrowsum = rowsums[j]
            rowminsum = float(sum(where(r1<r2, r1,r2)))
            if (irowsum == 0.0 and jrowsum == 0.0):
                cur_d = 0.0 # => two rows of zeros
            elif (irowsum == 0.0 or jrowsum == 0.0):
                cur_d = 1.0 # one row zeros, one not all zeros
            else:
                cur_d = 1.0 - (((rowminsum/irowsum) + (rowminsum/jrowsum))/2.0)
            dists[i][j] = dists[j][i] = cur_d
    return dists

def dist_manhattan(datamtx, strict=True):
    """ returns manhattan (city block) distance between rows

    dist(a,b) = sum on i( abs(a_i - b_i) )
    negative values ok (but not tested)
    
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * if strict==True, raises ValueError if any of the input data is
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a 2d matrix.  
    If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    for i in range(numrows):
        r1 = datamtx[i] # cache here
        for j in range(i):
            dists[i,j] = dists[j,i] = sum(abs(r1 - datamtx[j]))
            
    return dists

def dist_abund_jaccard(datamtx, strict=True):
    """Calculate abundance-based Jaccard distance between rows

    The abundance-based Jaccard index is defined in Chao et. al.,
    Ecology Lett. 8, 148 (2005), eq. 5:

    J_abd = UV / (U + V - UV), where

    U = sum of relative abundances of shared species in a
    V = sum of relative abundances of shared species in b

    The Chao-Jaccard distance is 1 - J_abd
    
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * an all zero row compared with a not all zero row returns a distance of 1
    * if strict==True, raises ValueError if any of the input data is negative,
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a matrix with nonnegative 
    entries.  If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    
    rowsums = datamtx.sum(axis=1, dtype='float')

    for i in range(numrows):
        row1 = datamtx[i]
        N1 = rowsums[i]

        for j in range(i):
            row2 = datamtx[j]
            N2 = rowsums[j]

            if N1 == 0.0 and N2 == 0.0:
                similarity = 1.0
            elif N1 == 0.0 or N2 == 0.0:
                similarity = 0.0
            else:
                shared = logical_and(row1, row2)
                u = sum(row1[shared]) / N1
                v = sum(row2[shared]) / N2
                # Verified by graphical inspection
                if u == 0.0 and v == 0.0:
                    similarity = 0.0
                else:
                    similarity = (u * v) / (u + v - (u * v))

            dists[i][j] = dists[j][i] = 1 - similarity
    return dists

def dist_morisita_horn(datamtx, strict=True):
    """ returns morisita-horn distance between rows

    dist(a,b) = 1 - 2*sum(a_i * b_i) /( (d_a + d_b)* N_a * N_b )

    see book: magurran 2004 pg 246
    
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * an all zero row compared with a not all zero row returns a distance of 1
    * if strict==True, raises ValueError if any of the input data is negative,
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a matrix with nonnegative 
    entries.  If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    
    rowsums = datamtx.sum(axis=1, dtype='float')
    row_ds = (datamtx**2).sum(axis=1, dtype='float') # these are d_a, etc
    
    for i in range(numrows):
        if row_ds[i] !=0.:
            row_ds[i] = row_ds[i] / rowsums[i]**2
    # this leaves row_ds zero if actually 0/0
    for i in range(numrows):
        row1 = datamtx[i]
        N1 = rowsums[i]
        d1 = row_ds[i]
        for j in range(i):
            row2 = datamtx[j]
            N2 = rowsums[j]
            d2 = row_ds[j]
            if N2 == 0.0 and N1==0.0:
                dist = 0.0
            elif N2 == 0.0 or N1==0.0:
                dist = 1.0
            else:
            # d's zero only if N's zero, and we already checked for that
                similarity = 2*sum(row1*row2)
                similarity = similarity / ( (d1 + d2) * N1 * N2 )
                dist = 1 - similarity
            dists[i][j] = dists[j][i] = dist
    return dists

def dist_pearson(datamtx, strict=True):
    """ Calculates pearson distance (1-r) between rows

    
    note that the literature sometimer refers to the pearson dissimilarity
    as (1 - r)/2 (e.g.: BC Blaxall et al. 2003: Differential Myocardial Gene
    Expression in the Development and Rescue of Murine Heart Failure)
    
    for pearson's r, see for example: Thirteen Ways to Look at the 
    Correlation Coefficient by J rodgers, 1988

    * distance varies between 0-2, inclusive.  
    * Flat rows (all elements itentical) will return a distance of 1 relative
    to any non-flat row, and a distance of zero to another flat row
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * an all zero row compared with a not all zero row returns a distance of 1
    * if strict==True, raises ValueError if any of the input data is
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a 2d matrix
    If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')

    rowmeans =  mean(datamtx, axis=1)
    rowstds =  std(datamtx, axis=1)

    dists = zeros((numrows,numrows),'d')
    n = float(numrows)

    for i in range(numrows):
        r1 = datamtx[i,:]
        r1m = rowmeans[i]
        r1dev = r1 - r1m
        for j in range(i):
            r2 = datamtx[j,:]
            r2m = rowmeans[j]
            r2dev = r2 - r2m
            top = sum(r1dev*r2dev)
            sum1 = sum(r1dev**2)
            sum2 = sum(r2dev**2)

            if (sum1 == 0.0 and sum2 == 0.0):
                r = 1.0
            elif (sum1 == 0.0 or sum2 == 0.0):
                r = 0.0
            else:
                bottom = sqrt(sum1 * sum2)
                r = top/bottom
                
            dists[i][j] = dists[j][i] = 1.0 - r
            
    return dists

def dist_soergel(datamtx, strict=True):
    """ Calculate soergel distance between rows of a matrix
    
    see for example Evaluation of Distance Metrics..., Fechner 2004
    dist(a,b) = sum on i( abs(a_i - b_i) ) / sum on i( max(a_i, b_i) )
    
    returns: a symmetric distance matrix, numrows X numrows
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
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    for i in range(numrows):
        r1 = datamtx[i,:]
        for j in range(i):
            r2 = datamtx[j,:]
            top = float(sum(abs(r1 - r2)))
            bot = float(sum(where(r1>r2, r1,r2)))
            if bot <= 0.0:
                cur_d = 0.0
            else:
                cur_d = top/bot
            dists[i][j] = dists[j][i] = cur_d
    return dists

def dist_spearman_approx(datamtx, strict=True):
    """ Calculate spearman rank distance (1-r) using an approximation formula
    
    considers only rank order of elements in a row, averaging ties 
    [19.2, 2.1, 0.03, 2.1] -> [3, 1.5, 0, 1.5]
    then performs dist(a,b) = 6 * sum(D^2) / (N*(N^2 - 1))
    where D is difference in rank of element i between row a and row b,
    N is the length of either row
    
    * formula fails for < 2 columns, returns a zeros matrix
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * if strict==True, raises ValueError if any of the input data is
    not finite, or if the input data is not a rank 2 array (a matrix), or if
    there are less than 2 colunms
    * if strict==False, assumes input data is a 2d matrix.
    If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
        if numcols < 2:
            raise ValueError("input matrix has < 2 colunms")
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    
    if numcols < 2:
        return dists # formula fails for < 2 elements per row
    
    for i in range(numrows):
        r1 = datamtx[i,:]
        rank1 = _rankdata(r1)
        for j in range(i):
            r2 = datamtx[j,:]
            rank2 = _rankdata(r2)
            rankdiff = rank1 - rank2
            dsqsum = sum((rankdiff)**2)
            dist = 6*dsqsum / float(numcols*(numcols**2-1))
            
            dists[i][j] = dists[j][i] = dist
    return dists

def dist_specprof(datamtx, strict=True):
    """returns a row-row species profile distance matrix
    
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * an all zero row compared with a not all zero row returns a distance of 1
    * if strict==True, raises ValueError if any of the input data is negative,
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a matrix with nonnegative 
    entries.  If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    for i in range(numrows):
        r1 = datamtx[i]
        r1sum = sum(r1)
        for j in range(i):
            r2 = datamtx[j]
            r2sum = sum(r2)
            if r1sum == 0.0 or r2sum == 0.0:
                if r1sum == 0.0 and r2sum == 0.0:
                    dist = 0.0
                else: dist = 1.0
            else:
                dist = norm((r1/r1sum) - (r2/r2sum))
            dists[i,j] = dists[j,i] = dist

    return dists

def binary_dist_otu_gain(otumtx):
    """ Calculates number of new OTUs observed in sample A wrt sample B
    
        This is an non-phylogenetic distance matrix analagous to unifrac_g. 
        The number of OTUs gained in each sample is computed with respect to
        each other sample.
    
    """
    result = []
    for i in otumtx:
        row = []
        for j in otumtx:
            gain = 0
            for i_val, j_val in zip(i, j):
                if i_val > 0 and j_val == 0:
                    gain += 1
            row.append(gain)
        result.append(row)
    return array(result)

def binary_dist_chisq(datamtx, strict=True):
    """Calculates binary chi-square dist between rows, returns dist matrix.

    converts input array to bool, then uses dist_chisq
    """
    datamtx = datamtx.astype(bool)
    datamtx = datamtx.astype(float)
    return dist_chisq(datamtx, strict=True)

def binary_dist_chord(datamtx, strict=True):
    """Calculates binary chord dist between rows, returns dist matrix.

    converts input array to bool, then uses dist_chisq
    for binary data, this is identical to a binary hellinger distance
    """
    datamtx = datamtx.astype(bool)
    datamtx = datamtx.astype(float)
    return dist_chord(datamtx, strict=True)

def binary_dist_sorensen_dice(datamtx, strict=True):
    """Calculates Sorensen-Dice distance btw rows, returning distance matrix.

    Note: Treats array as bool. This distance = 1 - dice's coincidence index
    see Measures of the Amount of Ecologic Association Between Species
    Author(s): Lee R. Dice, 1945
    The 'o' in sorensen should be a non-ascii char, but isn't here for ease
    of use
    this is identical to a binary bray-curtis distance, as well as the
    binary whittaker distance metric.

    a = num 1's in a
    b = num 1's in b
    c = num that are 1's in both a and b
    Dice dist = 1 - (2*c)/(a + b).
    also known as whittaker:
    whittaker = (a + b - c)/( 0.5*(a+b) ) - 1
    
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * negative input values are not allowed - will return nonsensical results 
    and/or throw errors
    """
    datamtx = datamtx.astype(bool)
    datamtx = datamtx.astype(float)
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    rowsums = datamtx.sum(axis=1)
    
    for i in range(numrows):
        row1 = datamtx[i]
        for j in range(i):
            row2 = datamtx[j]
            bottom = float(rowsums[i] + rowsums[j])
            cur_d = 0.0
            if bottom:
                cur_d = 1-(2*logical_and(row1,row2).sum()/bottom)
            dists[i][j] = dists[j][i] = cur_d
    return dists

def binary_dist_euclidean(datamtx, strict=True):
    """Calculates binary euclidean distance between rows, returns dist matrix.

    converts input array to bool, then uses dist_euclidean
    """
    datamtx = datamtx.astype(bool)
    datamtx = datamtx.astype(float)
    return dist_euclidean(datamtx, strict=True)

def binary_dist_hamming(datamtx, strict=True):
    """Calculates hamming distance btw rows, returning distance matrix.

    Note: Treats array as bool. 
    see for example wikipedia hamming_distance, 20 jan 2008

    hamming is identical to binary manhattan distance

    Binary hamming: 
    a = num 1's in a
    b = num 1's in b
    c = num that are 1's in both a and b
    hamm = a + b - 2c
    
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
    datamtx = datamtx.astype(bool)
    datamtx = datamtx.astype(float)
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    rowsums = datamtx.sum(axis=1)
    
    for i in range(numrows):
        first = datamtx[i]
        a = rowsums[i]
        for j in range(i):
            second = datamtx[j]
            b = rowsums[j]
            c = float(logical_and(first, second).sum())
            dist = a + b - (2.0*c)
            dists[i][j] = dists[j][i] = dist
    return dists
    
def binary_dist_jaccard(datamtx, strict=True):
    """Calculates jaccard distance between rows, returns distance matrix.

    converts matrix to boolean.  jaccard dist = 1 - jaccard index
    
    see for example: wikipedia jaccard index (20 jan 2009)
    this is identical to a binary version of the soergel distance

    Binary jaccard:
    a = num 1's in a
    b = num 1's in b
    c = num that are 1's in both a and b
    jaccard = 1 - (c/(a+b-c))
    
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
    datamtx = datamtx.astype(bool)
    datamtx = datamtx.astype(float)
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    
    rowsums = datamtx.sum(axis=1)
    for i in range(numrows):
        first = datamtx[i]
        a = rowsums[i]
        for j in range(i):
            second = datamtx[j]
            b = rowsums[j]
            c = float(logical_and(first, second).sum())
            if a==0.0 and b==0.0:
                dist = 0.0
            else:
                dist = 1.0 - (c/(a+b-c))
            dists[i][j] = dists[j][i] = dist
    return dists

def binary_dist_lennon(datamtx, strict=True):
    """Calculates lennon distance between rows, returns distance matrix.

    converts matrix to boolean.  jaccard dist = 1 - lennon similarity
    lennon's similarity is a modification of simpson's index
    see Jack J.  Lennon, The geographical structure of British bird 
    distributions: diversity, spatial turnover and scale

    Binary lennon:
    a = num 1's in a
    b = num 1's in b
    c = num that are 1's in both a and b
    lennon = 1 - (c/(c + min(a-c,b-c)))
    
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
    datamtx = datamtx.astype(bool)
    datamtx = datamtx.astype(float)
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    
    rowsums = datamtx.sum(axis=1)
    for i in range(numrows):
        first = datamtx[i]
        a = rowsums[i]
        for j in range(i):
            second = datamtx[j]
            b = rowsums[j]
            c = float(logical_and(first, second).sum())
            if a==0.0 and b==0.0:
                dist = 0.0
            elif c==0.0:
                dist = 1.0
            else:
                dist = 1.0 - (c/(c + min([a-c,b-c])))
            dists[i][j] = dists[j][i] = dist
    return dists

def binary_dist_ochiai(datamtx, strict=True):
    """Calculates ochiai distance btw rows, returning distance matrix.

    Note: Treats array as bool. 
    see for example:
    On the Mathematical Significance of the Similarity Index of Ochiai...
    Bolton, 1991

    a = num 1's in a
    b = num 1's in b
    c = num that are 1's in both a and b
    ochiai = 1 - (c/sqrt(a*b))
    
    * comparisons are between rows (samples)
    * input: 2D numpy array.  Limited support for non-2D arrays if 
    strict==False
    * output: numpy 2D array float ('d') type.  shape (inputrows, inputrows)
    for sane input data
    * two rows of all zeros returns 0 distance between them
    * an all zero row compared with a not all zero row returns a distance of 1
    * if strict==True, raises ValueError if any of the input data is negative,
    not finite, or if the input data is not a rank 2 array (a matrix).
    * if strict==False, assumes input data is a matrix with nonnegative 
    entries.  If rank of input data is < 2, returns an empty 2d array (shape:
    (0, 0) ).  If 0 rows or 0 colunms, also returns an empty 2d array.
    """
    datamtx = datamtx.astype(bool)
    datamtx = datamtx.astype(float)
    if strict:
        if not all(isfinite(datamtx)):
            raise ValueError("non finite number in input matrix")
        if any(datamtx<0.0):
            raise ValueError("negative value in input matrix")
        if rank(datamtx) != 2:
            raise ValueError("input matrix not 2D")
        numrows, numcols = shape(datamtx)
    else:
        try:
            numrows, numcols = shape(datamtx)
        except ValueError:
            return zeros((0,0),'d')

    if numrows == 0 or numcols == 0:
        return zeros((0,0),'d')
    dists = zeros((numrows,numrows),'d')
    rowsums = datamtx.sum(axis=1)
    
    for i in range(numrows):
        first = datamtx[i]
        a = rowsums[i]
        for j in range(i):
            second = datamtx[j]
            b = rowsums[j]
            c = float(logical_and(first, second).sum())
            if a==0.0 and b==0.0:
                dist = 0.0
            elif a==0.0 or b==0.0:
                dist = 1.0
            else:
                dist = 1.0 - (c/sqrt(a*b))
            dists[i][j] = dists[j][i] = dist
    return dists
    
def binary_dist_pearson(datamtx, strict=True):
    """Calculates binary pearson distance between rows, returns distance matrix

    converts input array to bool, then uses dist_pearson
    """

    datamtx = datamtx.astype(bool)
    datamtx = datamtx.astype(float)
    return dist_pearson(datamtx, strict=True)


if __name__ == "__main__":
    """ just a test run"""
    matrix1 = array(    [   [10,8,4,1],
                            [8,6,2,1],
                            [0,0,0,0],
                            [0,0,1,0],
                            [1,1,0,1],
                            [1,0,8,10],
                            [0,0,0,0],
                            [8,6,2,1],
                                    ])
    
    res = dist_euclidean(matrix1)
    print "euclidean distance result: \n"
    print res

