#!/usr/bin/env python
# File created on 30 Nov 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from qiime.util import MetadataMap

def merge_mapping_files(mapping_files, no_data_value='no_data'):
    """ Merge list of mapping files into a single mapping file 
        mapping_files: open file objects containing mapping data
        no_data_value: value to be used in cases where there is no
        mapping field associated with a sample ID (default: 'no_data')
    """
    metadata_maps = map(MetadataMap.parseMetadataMap, mapping_files)
    merged = sum(metadata_maps[1:], metadata_maps[0])
    merged.no_data_value = no_data_value
    return merged
