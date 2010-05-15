#!/usr/bin/env python
# File created on 21 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from random import shuffle
from numpy import array, mean
from cogent.cluster.procrustes import procrustes
from cogent.util.dict2d import Dict2D
from qiime.parse import parse_coords
from qiime.format import format_coords

def shuffle_full_matrix(m):
    """  """
    shape = m.shape
    data = m.flatten()
    shuffle(data)
    return array(data).reshape(shape)
    
def shuffle_within_rows(m):
    result = []
    for k in m:
        l = list(k)
        shuffle(l)
        result.append(l)
    return array(result)
    
def shuffle_col_order(m):
    order = range(m.shape[1])
    shuffle(order)
    result = []
    current_row = []
    for k in m:
        for o in order:
            current_row.append(k[o])
        result.append(current_row)
        current_row = []
    return array(result)


def map_sample_ids(sample_ids,sample_id_map):
    """ Map sample ids to new values in sample_id_map """
    try:
        result = [sample_id_map[sample_id] for sample_id in sample_ids]
    except KeyError:
        raise KeyError, "Unknown sample ID: %s" % sample_id
        
    return result

def reorder_coords(coords,sample_ids,order):
    """ Arrange the rows in coords to correspond to order 
    
        Note: order is the master list here -- if a sample id is
        not included in order, that coord will not be included in
        the results. All sample ids in order however must be in
        sample_ids
    
    """
    try:
        result =  array([coords[sample_ids.index(sample_id)]\
                    for sample_id in order])
    except ValueError:
        raise ValueError, 'Unknown sample ID: %s' % sample_id
    return result

def filter_coords_matrix(coords,dimensions_to_keep):
    return coords[:,:dimensions_to_keep]

def pad_coords_matrix(coords,dimensions_to_add):
    if dimensions_to_add < 0:
        raise ValueError, 'Dimensions to add must be >= 0.'
    elif dimensions_to_add == 0:
        return coords
    else:
        result = []
        # find a more efficent way to do this (don't have
        # internet right now)
        return array([list(c) + [0.] * dimensions_to_add
                      for c in coords])

def pad_coords_matrices(coords1,coords2):
    # Determine how many dimensions are in each vector, and the difference
    coords1_len = coords1.shape[1]
    coords2_len = coords2.shape[1]
    coords_len_diff = coords1_len - coords2_len
    # If the first vector is shorter, pad it with zeros
    if coords_len_diff < 0:
        coords1 = pad_coords_matrix(coords1,-1 * coords_len_diff)
    # If the second vector is shorter, pad it with zeros        
    elif coords_len_diff > 0:
        coords2 = pad_coords_matrix(coords2,coords_len_diff)
    else:
        pass
    return coords1, coords2
    
def get_mean_percent_variation(v1,v2):
    return [mean([p1,p2]) for p1,p2 in zip(v1,v2)]

def get_mean_eigenvalues(v1,v2):
    return [mean([p1,p2]) for p1,p2 in zip(v1,v2)]

def get_procrustes_results(coords_f1,coords_f2,sample_id_map=None,\
    randomize=None,max_dimensions=None,\
    get_eigenvalues=get_mean_eigenvalues,\
    get_percent_variation_explained=get_mean_percent_variation):
    """ """
    # If this is a random trial, apply the shuffling function passed as 
    # randomize()
    # Parse the PCoA files
    sample_ids1, coords1, eigvals1, pct_var1 = parse_coords(coords_f1)
    sample_ids2, coords2, eigvals2, pct_var2 = parse_coords(coords_f2)
    
    if sample_id_map:
        sample_ids1 = map_sample_ids(sample_ids1,sample_id_map)
        sample_ids2 = map_sample_ids(sample_ids2,sample_id_map)
    # rearrange the order of coords in coords2 to correspond to 
    # the order of coords in coords1 
    order = list(set(sample_ids1) & set(sample_ids2)) 
    coords1 = reorder_coords(coords1,sample_ids1,order)
    coords2 = reorder_coords(coords2,sample_ids2,order)
    
    if randomize:
        coords2 = randomize(coords2)
        
    coords1, coords2 = pad_coords_matrices(coords1,coords2)
    if max_dimensions:
        coords1 = filter_coords_matrix(coords1,max_dimensions)
        coords2 = filter_coords_matrix(coords2,max_dimensions)
        pct_var1 = pct_var1[:max_dimensions]
        pct_var2 = pct_var2[:max_dimensions]
    
    # Run the Procrustes analysis
    transformed_coords_m1, transformed_coords_m2, m_squared =\
     procrustes(coords1,coords2)
    
    eigvals = get_eigenvalues(eigvals1, eigvals2)
    pct_var = get_percent_variation_explained(pct_var1,pct_var2)
    
    transformed_coords1 = format_coords(coord_header=order,\
                                        coords=transformed_coords_m1,\
                                        eigvals=eigvals,\
                                        pct_var=pct_var)
    transformed_coords2 = format_coords(coord_header=order,\
                                        coords=transformed_coords_m2,\
                                        eigvals=eigvals,\
                                        pct_var=pct_var)
    
    # Return the results
    return transformed_coords1, transformed_coords2, m_squared

def procrustes_monte_carlo(coords_f1,\
                           coords_f2,\
                           trials=1000,\
                           max_dimensions=None,\
                           shuffle_f=shuffle_within_rows,\
                           sample_id_map=None):
    """ Run procrustes analysis with random trials
    
        This analysis could be made more efficient, as the current version 
        just calls get_procrustes_results() random_trials times, which involves
        re-parsing each time, etc. It is very fast though, so that's low 
        priority.
    
    """
    # Get the M^2 for the actual data
    actual_m_squared = get_procrustes_results(\
     coords_f1,\
     coords_f2,\
     sample_id_map=sample_id_map,\
     randomize=None,\
     max_dimensions=max_dimensions)[2]
        
    # Get the M^2 for the random trials, and count how many
    # are lower than or equal to the actual M^2
    trial_m_squareds = []
    count_better = 0
    for i in range(trials):
        trial_m_squared = get_procrustes_results(\
             coords_f1,\
             coords_f2,\
             sample_id_map=sample_id_map,\
             randomize=shuffle_f,\
             max_dimensions=max_dimensions)[2]
        trial_m_squareds.append(trial_m_squared)
        
        if trial_m_squared <= actual_m_squared:
            count_better += 1
    
    return actual_m_squared, trial_m_squareds, count_better, count_better/trials
