#!/usr/bin/env python
# File created on 13 Aug 2013
from __future__ import division

__author__ = "Luke Ursell"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Will Van Treuren, Luke Ursell"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Luke Ursell"
__email__ = "lkursell@gmail.com"
__status__ = "Development"

from biom.parse import parse_biom_table
from qiime.parse import parse_mapping_file_to_dict
from numpy import array
from qiime.pycogent_backports.test import ANOVA_one_way, correlation_test
from cogent.maths.stats.util import Numbers

from collections import defaultdict

"""To Do:

"""

otu_table_fp = "/Users/lukeursell/Desktop/otu_table.biom"
mapping_file_fp = '/Users/lukeursell/Desktop/Fasting_Map.txt'

# def get_biom_data(otu_table_fp):
#     """parse out biom table, get important info 
#     """
#     # get biom table data in array format    
#     bt = parse_biom_table(open(otu_table_fp))
#     bt_data = array([bt.observationData(i) for i in bt.ObservationIds])

#     #get list of OTU ids for writing out results
#     otu_ids = [i for i in bt.ObservationIds]
    
#     return bt_data, otu_ids

# might move this to script level...
# def parse_mapping(mapping_file_fp):
#     """Parses mapping file

#     Input: mapping file filepath

#     Output:
#     map_data is a list of all metadata:
#     ['PC.636',
#       'ACGGTGAGTGTC',
#       'YATGCTGCCTCCCGTAGGAGT',
#       'Fast',
#       '',
#       '20080116',
#       'Fasting_mouse_I.D._636']

#     map_headers is a list of headers:
#       ['SampleID',
#      'BarcodeSequence',
#      'LinkerPrimerSequence',
#      'Treatment',
#      'test_col',
#      'DOB',
#      'Description']
#     """
#     map_data, map_headers, _ = parse_mapping_file(mapping_file_fp)

#     return map_data, map_headers

## These two functions could probably be combined somehow...

# def get_category_info(map_data, map_headers, category):
#     """Create a dict of {SampleID: category_value} and continuous category values 

#     When the category of interest is catageorical, the dicionary produced
#     contains the SampleId and which categorical group it belongs to.
#     When the category contains continuous data, the category_values is produced
#     and contains a list of all values
    
#     Input: map_data, map_headers (from parse_mapping), and category of interest
#     e.g. cat_info, category_values = get_category_info(map_data, map_headers, 'Treatment')

#     Output: 
#     cat_info: dictionary relating SampleIds
#     {'PC.354': 'Control',
#      'PC.355': 'Control',
#      'PC.356': 'Control',
#      'PC.481': 'Control',
#      'PC.593': 'Control',
#      'PC.607': 'Fast',
#      'PC.634': 'Fast',
#      'PC.635': 'Fast',
#      'PC.636': 'Fast'}

#      category_values: ['10', '20', '30', '40', '50', '60', '70', '80', '90']
#     """
#     cat_info = {}
#     category_values = []

#     # find index in mapping data of category of interest
#     category_index = map_headers.index(category)
    
#     # walk through each SampleID individually
#     for line in map_data:
#         sample_id = line[0]
#         category_val = line[category_index]
        
#         # if the category value is blank in the mapping file, ignore SampleID
#         if category_val != "":
#             cat_info[sample_id] = category_val
#             if category_val not in category_values:
#                 category_values.append(category_val)
#         elif category_val == "":
#             print "Sample %s contained an empty field for the category %s \
#                 and will be ignored" % (sample_id,category)
#             pass

#     return cat_info, category_values

# this function could probably be incorporated into get_category_info
# feel free to combine them if you see an easy way



# def build_category_sampleid_dict(category_result):
#     """Take output of get_category_info and build a {category: [list of ids]}
#     like: 
#     {'Control': ['PC.355', 'PC.593', 'PC.356', 'PC.481', 'PC.354'],
#      'Fast': ['PC.636', 'PC.607', 'PC.634', 'PC.635']}
#     """
#     cat_to_ids = {}
#     for k, v in category_result.iteritems():
#         cat_to_ids.setdefault(v, []).append(k)
#     return cat_to_ids


# this would be run at the beginning of the script



# def get_sampleid_indices(cat_to_ids, otu_table_fp):
#     """Take in put of cat_to_ids, and replace the SampleIDs with their
#     indexed position in the otu table:

#     Output:
#     {'Control': [6, 5, 2, 3, 4], 'Fast': [0, 7, 8, 1]}
#     """
#     #Parse otu table and create list of SampleIds found in table
#     otu_table = parse_biom_table(open(otu_table_fp))
#     otu_table_ids = [i for i in otu_table.SampleIds]

#     cat_to_sampleid_index = {}
    
#     for cat,ids in cat_to_ids.iteritems():
#         # check if category is already added
#         if cat not in cat_to_sampleid_index:
#             cat_to_sampleid_index[cat] = []
            
#             # append SampleID indexed position in otu table to dict
#             for id in ids:
                
#                 #check to see if SampleID is in the mapping file but not the
#                 # OTU table
#                 if id not in otu_table_ids:
#                     print "Sample %s is not found in the otu table and will \
#                         be ignored" % id
                
#                 #if SampleID is in both mapping file and OTU table, get index
#                 else: 
#                     id in otu_table_ids        
#                     index = otu_table.getSampleIndex(id)
#                     cat_to_sampleid_index[cat].append(index)
#         else:
#             continue

#     return cat_to_sampleid_index

# def get_category_arrays(cat_to_sampleid_index, bt_data):
#     """Take bt_data array, and pull out indexed positions for each category of
#     interest.

#     Returns a list of arrays:
#     [[array([[ 0.,  0.,  0.,  0.,  0.],
#        [ 0.,  0.,  0.,  0.,  0.],
#        ...
#        [ 0.,  0.,  0.,  0.,  0.]])],
#     [array([[ 1.,  0.,  0.,  0.],
#        [ 0.,  0.,  0.,  1.],
#        ..., 
#        [ 0.,  0.,  0.,  0.],
#        [ 1.,  0.,  0.,  0.]])]]
#     """
#     listed_category_data = []
#     for cat in cat_to_sampleid_index.keys():
#         listed_category_data.append(bt_data[:,cat_to_sampleid_index[cat]])

#     return listed_category_data


def sync_biom_and_mf(pmf, bt):
    """Reduce mapping file dict and biom table to shared samples."""
    mf_samples = set(pmf.keys())
    bt_samples = set(bt.SampleIds)
    if mf_samples == bt_samples:
        # agreement, can continue without fear of breaking code
        pass
    else: 
        shared_samples = mf_samples.intersection(bt_samples)
        # check that we shared something
        assert len(shared_samples)!=0, \
            "sync_biom_and_mf: No shared samples, no point in continuing."
        # tell the user which samples we are excluding
        print "The following samples were not shared, and will not be "+\
            "considered in the analysis:\n" + \
            ', '.join(shared_samples-(mf_samples.union(bt_samples)))
        # remove samples that were in the mapping file but not biom file
        [pmf.pop(i) for i in pmf if i not in shared_samples]
        # remove samples in the biom table that were not in the mapping file
        def _f(sv, sid, smd):
            if sid in shared_samples:
                return True
            else:
                return False
        bt = bt.filterSamples(_f)
    return pmf, bt

def get_sample_cats(pmf, category):
    """Create {SampleID:category_value} for samples in parsed mf dict."""
    return {k:pmf[k][category] for k in pmf.keys() if pmf[k][category] != ""}

def get_cat_sample_groups(sam_cats):
    """Create {category_value:[samples_with_that_value} dict."""
    cat_sam_groups = {group:[] for group in set(sam_cats.values())}
    [cat_sam_groups[v].append(k) for k,v in sam_cats.items()]
    return cat_sam_groups

def get_sample_indices(cat_sam_groups, bt):
    """Create {category_value:index_of_sample_with_that_value} dict."""
    return {k:[bt.SampleIds.index(i) for i in v] for k,v in cat_sam_groups.items}

def row_generator(bt, cat_sam_groups):
    """Produce a generator that can feed lists of arrays to any test."""
    data = array([bt.observationData(i) for i in bt.ObservationIds])
    return ([row[cat_sam_groups[k]] for k in cat_sam_groups] for row in data)

def run_ANOVA(listed_category_data):
    """Compute the ANOVA for inputed lists of data

    Take in the list of array data produced by get_category_arrays.
    Runs ANOVA_one_way on each OTU.

    Returns a list of group means and the propability from ANOVA

        [([0.2, 0.0], 0.40708382206558885),
        ([0.0, 0.25], 0.29235199244023835)]
    """
    result = []
    counter = 0

    while counter < len(listed_category_data[0]):
        a = []
        
        # Create a list of lists where contents are in Numbers format
        a.append(Numbers(listed_category_data[0][counter].tolist()))
        a.append(Numbers(listed_category_data[1][counter].tolist()))
        
        # Compute ANOVA
        dfn, dfd, F, between_MS, within_MS, group_means, prob = ANOVA_one_way(a)
        
        # Append ANOVA output to result
        output = group_means, prob
        result.append(output)
        counter += 1

    return result

def run_correlation(bt_data, category_values, corr_method, tails=None, 
                    permutations=999, confidence_level=0.95):
    """Compute specified correlation between x and y data

    Take in bt_data (an array of the biom table observations), a list of the
    category_values (obtained from the mapping file), and the desired
    correlation method (either 'spearman' or 'pearson')

    Returns list of results, each entry containing corr_coeff, parametric_p_val,
    and nonparametric_p_val
    """
    result = []
    
    #Convert continuous variables from the mapping file to integers
    cat_vals = [int(i) for i in category_values]

    for otu in bt_data:
        # compute correlation
        (corr_coeff, parametric_p_val, permuted_corr_coeffs, \
        nonparametric_p_val, (ci_low, ci_high)) = \
        correlation_test(otu, cat_vals, corr_method, tails, permutations, confidence_level)

        # grab selected correlation info and append to results
        output = corr_coeff, parametric_p_val, nonparametric_p_val
        result.append(output)

    return result







