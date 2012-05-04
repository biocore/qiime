#!/usr/bin/env python
# File created on 02 May 2012
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "William Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"
 
from numpy.random import shuffle
from numpy import std, mean, array
from copy import deepcopy


"""NOTE: calculates NRI/NTI according to the formula used by Phylocom 4.2,3.41 
rather than Webb 2002 or Webb 2000. Uses a 'null model 2' -- chooses the 
random distance matrices without replacment from the full available pool of
distance matrices. See Phylocom manual for details."""

# NRI

def nri(dist_mat, all_ids, group_ids, iters=1000):
    """calculates nri (net relatedness index) of an otu or sample grouping
    inputs:
        datamtx - 2d array that has distance between otus or samples
        all_ids - list of strs/floats/ints that corresponds to the rows or cols 
            of the datamtx
        group_ids - list of strs/floats/ints that corresponds to the rows or
            cols of the datamtx that you want to group and compare against the 
            entire matrix.
        iters - int, number of random draws used to calculate the distance to
            compare the group_ids distance against.
    NOTE: calculates NRI according to the formula used by Phylocom 4.2,3.41 
    rather than Webb 2002 or Webb 2000. Uses a 'null model 2' -- chooses the 
    random distance matrices without replacment from the full available pool of
    distance matrices. See Phylocom manual for details. 
    """
    if len(group_ids) < 2:
        raise ValueError('there is no standard deviation of distances with ' +\
            'less than 2 taxa in taxa ids')
    if dist_mat.sum() == 0.0:
        raise ValueError('the input dist_mat has no data. its sum is 0.')

    mn_x_obs = \
        mpd(take_distmat_data(dist_mat, all_ids, group_ids))
    mn_x_n, sd_x_n = \
        mpd_mean_sd(dist_mat, all_ids, len(group_ids), iters)
    # for debugging, check against phylocom
    # print mn_x_obs, mn_x_n, sd_x_n
    
    return -1.0*((mn_x_obs-mn_x_n)/sd_x_n)

def take_random_ids(all_ids, num_to_take):
    """takes num_to_take ids from shuffled list of all_ids"""
    if len(all_ids) < num_to_take:
        raise ValueError('trying to take too many ids from all ids')
    l = deepcopy(all_ids)
    shuffle(l)
    return l[:num_to_take]

def mpd(datamtx):
    """calculates mpd (mean phylogenetic distance)""" 
    dists = []
    rows,cols = datamtx.shape
    for r in range(rows):
        for c in range(r+1,rows): # avoid d(i,i). d(i,i)=0 if datamtx calculated with distance metric
            dists.append(datamtx[r][c])
    m = mean(dists)
    return float(m)

def mpd_mean_sd(dist_mat, all_ids, num_to_take, iters):
    """calculates mean and std for random dist mats based on iter randomizations
    """
    means = []
    for i in range(iters):
        ids_to_take = take_random_ids(all_ids, num_to_take)
        r_dist_mat = take_distmat_data(dist_mat, all_ids, ids_to_take)
        i_mean = mpd(r_dist_mat)
        means.append(i_mean)
    # calculate the standard deviation of the means of the distances rather than
    # the standard deviation of the distances themselves. The forumula from
    # Webb 2002 seems to calculate the standard deviation of the distances
    # but based on tests of Phylocom and the Phylocom manual, Phylocom 
    # computes the standard deviations of the means of the iters number of 
    # random  distance matrices.  
    return float(mean(means)), float(std(means))

def take_distmat_data(dist_mat, all_ids, ids_to_keep):
    """takes data from the dist_mat based on index of ids_to_keep in all_ids"""
    # ids could be given in any order, must sort  to prevent choosing
    # wrong values at later steps. all_ids must be given in order of the
    # dist_mat array they correspond to
    indices = sorted([all_ids.index(i) for i in ids_to_keep])
    new_dist_mat = dist_mat.take(indices,axis=0) #remove rows
    new_dist_mat = new_dist_mat.take(indices,axis=1) #remove cols
    return new_dist_mat


# NTI

def nti(datamtx, all_ids, group_ids, iters=1000):
    """calculates nti (nearest taxon index) of an otu or sample grouping
    inputs:
        datamtx - 2d array that has distance between otus or samples
        all_ids - list of strs/floats/ints that corresponds to the rows or cols 
            of the datamtx
        group_ids - list of strs/floats/ints that corresponds to the rows or
            cols of the datamtx that you want to group and compare against the 
            entire matrix.
        iters - int, number of random draws used to calculate the distance to
            compare the group_ids distance against.
    NOTE: calculates NTI according to the formula used by Phylocom 4.2,3.41 
    rather than Webb 2002 or Webb 2000. Uses a 'null model 2' -- chooses the 
    random distance matrices without replacment from the full available pool of
    distance matrices. See Phylocom manual for details. 
    """
    if len(group_ids) < 2:
        raise ValueError('there is no standard deviation of distances with ' +\
            'less than 2 taxa in the taxa ids')
    if datamtx.sum() == 0.0:
        raise ValueError('the input dist_mat has no data. its sum is 0.')

    mn_y_obs = \
        mntd(take_distmat_data(datamtx, all_ids, group_ids))
    mn_y_n, sd_y_n = \
        mntd_mean_sd(datamtx, all_ids, len(group_ids), iters)
    # for debugging, to check against phylocom
    # print mn_y_obs, mn_y_n, sd_y_n
    
    return -1.*((mn_y_obs-mn_y_n)/sd_y_n)

def mntd(datamtx):
    """calculates mntd (mena nearest taxon distance) for a datamtx"""
    # the MNTD (mean nearest taxon distance) requires choosing only the min
    # value in each row, excluding d(i,i) vals which are dists to self (0.0)
    m = 0.0
    for i,row in enumerate(datamtx): 
        m += min(row.tolist()[:i]+row.tolist()[i+1:]) # exclude d(i,i). d(i,i)=0 bc datamatx calculated with metric   
    return m/float(i+1) #i+1 dists 


def mntd_mean_sd(dist_mat, all_ids, num_to_take, iters):
    """calculates mean and std for random dist mats based on iter randomizations
    num_to_take is the number of random rows/cols to take for a randomization"""
    means = []
    for i in range(iters):
        ids_to_take = take_random_ids(all_ids, num_to_take)
        r_dist_mat = take_distmat_data(dist_mat, all_ids, ids_to_take)
        i_mean = mntd(r_dist_mat)
        means.append(i_mean)
    return float(mean(means)), float(std(means))

