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

def get_biom_data():
	bt = parse_biom_table(open("/Users/lukeursell/Desktop/otu_table.biom"))
	bt_data = array([bt.observationData(i) for i in bt.ObservationIds])
	print bt_data

def parse_mapping():
	"""Parse mapping file
	"""
	map_data, map_headers, map_comments = \
		parse_mapping_file('/Users/lukeursell/Desktop/Fasting_Map.txt')

	return map_data, map_headers

def get_category_info(map_data, map_headers, category):
	"""Create a dict of {SampleID: category_value} and a list of all possible
	category_values (to be used when the category contains continuous data)
	"""
	result = {}
	category_values = []

	category_index = map_headers.index(category)
	for line in map_data:
		sample_id = line[0]
		category_val = line[category_index]
		
		# if the category value is blank in the mapping file, ignore SampleID
		if category_val != "":
			result[sample_id] = category_val
			if category_val not in category_values:
				category_values.append(category_val)
		elif category_val == "":
			continue

	return result, category_values

def build_category_sampleid_dict(category_result):
	cat_to_ids = {}
	for i in category_result.itervalues():
		if 










