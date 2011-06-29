#!/usr/bin/env python
# File created on 19 May 2011
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren","Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "William Van Treuren"
__email__ = "vantreur@colorado.edu"
__status__ = "Release"

from numpy import array, isnan
from qiime.parse import parse_mapping_file_to_dict, parse_rarefaction
from cogent.maths.stats.test import t_two_sample

def make_value_pairs_from_category(mapping_data, category):
    """creates all pairs of unique category values from mapping data
    
    input:
        mapping_data - a nested dictionary that maps SampleIds to
        descriptor categories (e.g. {'Id1': {'Weight':'Fat'}}
    
        category - a string which specifies a category found in the 
        mapping data (e.g. 'Obese')
    
    output:
        unique_pairs - a list of unique pairs of the values specified by
        category
        (e.g. [('Obese','Fat'),('Obese','notFat'),('Fat','notFat')]
    
    """
    
    #gets the keys of the mapping file dictionary. corresponds to the 
    #names of the individuals in the mapping data file
    
    keys = mapping_data.keys()
       
    categories = []
    for key in keys:
        try:
            categories.append(mapping_data[key][category])
        except KeyError:
            raise ValueError(('the specified category ({0}) was '+\
                    'not found in the mapping file.').format(category))
    #strip duplicate values from this list
    
    unique_vals = []
    for val in categories:
        if unique_vals.count(val)==0:
            unique_vals.append(val)
            
    #create and populate list of unique pairs of category values
    
    unique_pairs = []
    for i in range(len(unique_vals)):
        for j in unique_vals[i+1:]:
            unique_pairs.append((unique_vals[i], j))
    
    return unique_pairs
            
def make_category_values_Id_dict(mapping_data, category):
    """makes dict lists all SampleIds that have given category value
    
    input:
        mapping_data - a nested dictionary that maps SampleID to
        descriptor categories (e.g. {'Id1': {'Weight':'Fat}}
    
        category - a string which specifies a category found in the 
        mapping data (e.g. 'Weight')
    
    output:
        cat_val_Ids - a dictionary with category values as keys
        and a list of SampleIds which have the specific category value
        as values 
        (e.g. {'Fat':['Id1','Id2'],'notFat':['Id3'],'Obese':['Id4']}).
    
    """
    
    keys = mapping_data.keys()
    
    # create and populate list of all the different values for the 
    # category that was specified 
    
    categories = []
    for key in keys:
        categories.append(mapping_data[key][category])
    
    #strip duplicate values from this list
    
    unique_vals = []
    for val in categories:
        if unique_vals.count(val)==0:
            unique_vals.append(val)
    
    #make a dictionary with keys that are the possible values of the 
    #category that was specified.
    
    cat_val_Ids = {}
    for val in unique_vals:
        cat_val_Ids[val] = []
    
    #populate the cat_val dict with Id's which have proper category 
    
    for key in keys:
        for val in cat_val_Ids.keys():
            if mapping_data[key][category] == val:
                cat_val_Ids[val].append(key)
    
    return cat_val_Ids
        
def map_category_value_pairs_to_Ids(value_pairs, cat_val_Ids):
    """maps category value pairs to Id's which have that category value
    
    input:
        value_pairs - a list of pairs of categories (e.g. 
        [('Obese','Fat'),('Fat','notFat'),('Obese','notFat')]
        
        cat_val_Ids - a dictionary with category values as keys
        and a list of SampleId's which have the specific category value 
        as values
        (e.g. {'Fat':['Id1','Id2'],'notFat':['Id3'],'Obese':['Id4']}).
    
    output:
        mapped_pairs - the list of value_pairs with the values replaced
        by the SampleIds which have the values specified in the pair e.g
        [(['Id4'],['Id1','Id2']),(['Id1','Id2'],\
        ['Id3]),(['Id4'],['Id3])] 
    """
    
    mapped_pairs = []
    
    for pair in value_pairs:
        mapped_pairs.append((cat_val_Ids[pair[0]],cat_val_Ids[pair[1]]))
    
    return mapped_pairs


def make_SampleIds_rarefaction_columns_dict(rarefaction_list):
    """maps SampleId to column in parsed rarefaction file output
    
    input:
        rarefaction_list - ouput of parse_rarefaction.py. a nested list
        of scores and SampleIds and other fields. 
    
    output:
        map_from_Id_to_col - a dict which has as keys the SampleIds, and
        as values the col they are in the in the parsed rarefaction list
    """
    
    map_from_Id_to_col = {}
    
    # the first 3 entries in the rarefaction_list are not SampleIDs
    # so we dump them
    
    Ids = rarefaction_list[0][3:]
    
    for Id in Ids:
        map_from_Id_to_col[Id] = Ids.index(Id)
    
    return map_from_Id_to_col
        
    
def extract_rarefaction_scores_at_depth(depth, rarefaction_list):
    """makes rarefaction matrix with row=iteration and col=SampleId
    
    input:
        depth - an integer which corresponds to the depth of the
        rarefaction. also called the "sequences per sample" in the 
        rarefaction file
        
        rarefaction_list - ouput of parse_rarefaction.py. a nested list
        of scores and SampleIds and other fields. 
    
    output:
        result - a matrix with rows=rarefaction scores at a given depth
        and iteration, and cols=SampleIds.
    
    """
    # make and populate an array that has as rows rarefaction values
    # at the same depth and iteration and as cols SampleIds.
    
    score_matrix = []
    
    # the 4th element of rarefaction_list is a list of scores for each
    # different SampleId
    
    for line in rarefaction_list[3]:
        if line[0] == depth:
            # the first two elements are just rarefaction depth and 
            # iteration, throw these away
            score_matrix.append(line[2:])
    
    # raise error if rarefaction depth not found in rarefaction file
    if score_matrix == []:
        raise ValueError(('Specified depth ({0}) was not found in '+\
                    'the rarefaction file.').format(depth))
                        
    
    score_matrix_elements = []
    
    for line in score_matrix:
        score_matrix_elements.append(line)
    
    result = array(score_matrix_elements) 
    
    # raise error if any rarefaction score at spec. depth is Nan
    if isnan(result).any():
        raise ValueError(('Specified depth ({0}) has NaNs for some '+\
                            'rarefaction scores.').format(depth))
    
    return result


def convert_SampleIds_to_rarefaction_mtx(chosen_SampleIds,score_matrix,\
                                        map_from_SampleIds_to_cols):
    """converts list of SampleIDs to score mtx from rarefaction file
    
    input:
        chosen_SampleIds - a list of SampleIds
        
        score_matrix - a matrix created by
        extract_rarefaction_scores_at_depth which represents the
        rarefaction scores for a given depth
        
        map_from_SampleIds_to_cols - a dict which maps a SampleId to 
        the column its scores are in in the score matrix
    
    output:
        reduced_scores_matrix - a matrix which is the input scores mtx
        with only the cols that correspond to the chosen_SampleIds
    
    """
    #create and populate a list that specifies the r_array columns which
    #correspond to the name_list
    
    cols=[]
    
    for Id in chosen_SampleIds:
        cols.append(map_from_SampleIds_to_cols[Id])
        
    # grab only the columns we need based on a passed list of names and
    # a dictionary to convert between those names and the proper cols
    reduced_scores_matrix = score_matrix.take(cols, axis=1)
    
    return reduced_scores_matrix
    
    

def compare_alpha_diversities(rarefaction_lines, mapping_lines, 
                              category, depth):
    """compares alpha diversities
    
    inputs:
        rarefaction_file - rarefaction file which gives scores for 
        various rarefactions and depths
        
        mapping_file - file that has ID's and categories that the ID's
        fall in
        
        category - the category to be compared, is a string
        
        depth - the depth of the rarefaction_file to use, is an integer
    
    outputs:
        results - a nested dictionary which specifies the category as
        the top level key, and as its value, dictionaries which give the
        results of the t_two_sample test for all unique pairs of values
        in the specified category
    
    """
     
    rarefaction_data = parse_rarefaction(rarefaction_lines)
    mapping_data = parse_mapping_file_to_dict(mapping_lines)[0]
    value_pairs = make_value_pairs_from_category(mapping_data, category)
    
    category_values_Ids = make_category_values_Id_dict(mapping_data, 
                                                       category)
    
    SampleId_pairs = map_category_value_pairs_to_Ids(value_pairs,
                                                    category_values_Ids)
    
    map_from_Id_to_col = make_SampleIds_rarefaction_columns_dict(
                                                       rarefaction_data)
    
    reduced_rarefaction_mtx = extract_rarefaction_scores_at_depth(depth,
                                                       rarefaction_data)
    
    results = {category:{}}
    
    for pair in range(len(SampleId_pairs)):
        i=(convert_SampleIds_to_rarefaction_mtx(SampleId_pairs[pair][0],
                           reduced_rarefaction_mtx, map_from_Id_to_col))
        
        j=(convert_SampleIds_to_rarefaction_mtx(SampleId_pairs[pair][1],
                           reduced_rarefaction_mtx, map_from_Id_to_col))
        
        results[category][(str(value_pairs[pair][0]),
                           str(value_pairs[pair][1]))] =\
                          t_two_sample(i,j)
    
    return results
    
