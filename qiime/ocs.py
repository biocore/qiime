#!/usr/bin/env python
# File created on 13 Aug 2013
from __future__ import division

__author__ = "Luke Ursell"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Luke Ursell"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Luke Ursell"
__email__ = "lkursell@gmail.com"
__status__ = "Development"

from biom.parse import parse_biom_table
from qiime.parse import parse_mapping_file
from numpy import array

"""To Do:

"""

otu_table_fp = "/Users/lukeursell/Desktop/otu_table.biom"

def get_biom_data(otu_table_fp):
    
    bt = parse_biom_table(open(otu_table_fp))
    bt_data = array([bt.observationData(i) for i in bt.ObservationIds])
    return bt_data

def parse_mapping():
    """Parse mapping file

    map_data is a list of all metadata:
    ['PC.636',
      'ACGGTGAGTGTC',
      'YATGCTGCCTCCCGTAGGAGT',
      'Fast',
      '',
      '20080116',
      'Fasting_mouse_I.D._636']

  map_headers is a list of headers:
      ['SampleID',
     'BarcodeSequence',
     'LinkerPrimerSequence',
     'Treatment',
     'test_col',
     'DOB',
     'Description']
    """
    map_data, map_headers, map_comments = \
        parse_mapping_file('/Users/lukeursell/Desktop/Fasting_Map.txt')

    return map_data, map_headers

## These two functions could probably be combined somehow...
def get_category_info(map_data, map_headers, category):
    """Create a dict of {SampleID: category_value} and a list of all possible
    category_values (to be used when the category contains continuous data)

    e.g. result, category_values = get_category_info(map_data, map_headers, 'Treatment')

    results:
    {'PC.354': 'Control',
     'PC.355': 'Control',
     'PC.356': 'Control',
     'PC.481': 'Control',
     'PC.593': 'Control',
     'PC.607': 'Fast',
     'PC.634': 'Fast',
     'PC.635': 'Fast',
     'PC.636': 'Fast'}

     category_values: ['Control', 'Fast']
    """
    result = {}
    category_values = []

    # find index in mapping data of category of interest
    category_index = map_headers.index(category)
    
    # walk through each SampleID individually
    for line in map_data:
        sample_id = line[0]
        category_val = line[category_index]
        
        # if the category value is blank in the mapping file, ignore SampleID
        if category_val != "":
            result[sample_id] = category_val
            if category_val not in category_values:
                category_values.append(category_val)
        elif category_val == "":
            print "Sample %s contained an empty field for the category %s \
                and will be ignored" % (sample_id,category)
            pass

    return result, category_values

# this function could probably be incorporated into get_category_info
# feel free to combine them if you see an easy way
def build_category_sampleid_dict(category_result):
    """Take output of get_category_info and build a {category: [list of ids]}
    like: 
    {'Control': ['PC.355', 'PC.593', 'PC.356', 'PC.481', 'PC.354'],
     'Fast': ['PC.636', 'PC.607', 'PC.634', 'PC.635']}
    """
    cat_to_ids = {}
    for k, v in category_result.iteritems():
        cat_to_ids.setdefault(v, []).append(k)
    return cat_to_ids

def get_sampleid_indices(cat_to_ids, otu_table_fp):
    """Take in put of cat_to_ids, and replace the SampleIDs with their
    indexed position in the otu table:

    Output:
    {'Control': [6, 5, 2, 3, 4], 'Fast': [0, 7, 8, 1]}
    """
    #Parse otu table and create list of SampleIds found in table
    otu_table = parse_biom_table(open(otu_table_fp))
    otu_table_ids = [i for i in otu_table.SampleIds]

    cat_to_sampleid_index = {}
    
    for cat,ids in cat_to_ids.iteritems():
        # check if category is already added
        if cat not in cat_to_sampleid_index:
            cat_to_sampleid_index[cat] = []
            
            # append SampleID indexed position in otu table to dict
            for id in ids:
                
                #check to see if SampleID is in the mapping file but not the
                # OTU table
                if id not in otu_table_ids:
                    print "Sample %s is not found in the otu table and will \
                        be ignored" % id
                
                #if SampleID is in both mapping file and OTU table, get index
                elif: id in otu_table_ids:        
                    index = otu_table.getSampleIndex(id)
                    cat_to_sampleid_index[cat].append(index)
        else:
            continue

    return cat_to_sampleid_index

def get_category_arrays(cat_to_sampleid_index, bt_data):
    """Take bt_data array, and pull out indexed positions for each category of
    interest.

    Returns a list of arrays:
    [[array([[ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       ...
       [ 0.,  0.,  0.,  0.,  0.]])],
 [array([[ 1.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  1.],
       ..., 
       [ 0.,  0.,  0.,  0.],
       [ 1.,  0.,  0.,  0.]])]]
    """
    listed_category_data = []
    for cat in cat_to_sampleid_index.keys():
        listed_category_data.append([bt_data[:,cat_to_sampleid_index[cat]]])

    return listed_category_data











