#!/usr/bin/env python
# File created on 02 Nov 2012
from __future__ import division

__author__ = "Yoshiki Vazquez-Baeza"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Yoshiki Vazquez-Baeza"
__email__ = "yoshiki89@gmail.com"
__status__ = "Development"


from math import floor
from numpy import array
from copy import deepcopy

def check_mapping_file_headers():
    pass


def add_alpha_diversity_values_to_mapping_file(alpha_metrics, alpha_sample_ids,\
                                            alpha_data, mapping_file_headers,\
                                            mapping_file_data, number_of_bins,\
                                            missing_value_name='N/A'):
    """ """

    new_mapping_file_data = deepcopy(mapping_file_data)
    new_mapping_file_headers = deepcopy(mapping_file_headers)

    for index, metric in enumerate(alpha_metrics):
        # get the diversity value for this metric
        data = [[row[index]] for row in alpha_data]
        metric_max = max(data)[0]
        metric_min = min(data)[0]

        # add headers for each metric
        new_mapping_file_headers.append('{0}_alpha'.format(metric))
        new_mapping_file_headers.append('{0}_normalized_alpha'.format(metric))
        new_mapping_file_headers.append('{0}_alpha_label'.format(metric))

        # get the normalized value of diversity and the tag for each value
        for value in data:
            norm_value = (value[0]-metric_min)/(metric_max-metric_min)
            value.append(norm_value)
            value.append(_get_level(norm_value, number_of_bins, 'bin'))

        # iterate using the mapping file instead of using the alpha diversity
        # data because more often that you will like you will have more samples
        # in the mapping file than in your downstream analysis data
        for row in new_mapping_file_data:
            try:
                data_index = alpha_sample_ids.index(row[0])

                # data fields should be strings
                row.extend(map(str, data[data_index]))
            except ValueError:
                row.extend([missing_value_name, missing_value_name,\
                    missing_value_name])

    return new_mapping_file_data, new_mapping_file_headers

def _get_level(value, levels, prefix=None):
    """ """
    assert value <= 1, "The value cannot be greater than 1"
    assert value >= 0, "The value cannot be less than 0"
    assert levels > 0 and type(levels) is int, "The number of levels must be "+\
        "an integer and the value must be greater than zero."

    factor = 1/levels
    value_level = int(floor(value/factor))+1

    # take care of assignments where the division equals 1
    if value_level > levels:
        value_level = levels

    if prefix != None:
        output = '{0}_{1}_of_{2}'.format(prefix, value_level, levels)
    else:
        output = value_level

    return output

def alpha_diversity_data_to_dict(metrics, sample_ids, data):
    """ """

    assert len(sample_ids) == len(data), "There has to be as many sample_ids"+\
        "as rows in the alpha diversity data, cannot continue."
    assert len(metrics) == len(data[0]), "There has to be as many metrics as"+\
        "columns in the alpha diversity data, cannot continue."

    # iterate over the metrics to create dictionaries with the sample ids
    output = {}
    for metric_index, metric in enumerate(metrics):
        _buffer = {}
        for sample_id_index, sample_id in enumerate(sample_ids):
            _buffer[sample_id] = data[sample_id_index][metric_index]
        output[metric] = _buffer

    return output

def alpha_diversity_dict_to_data(dictionary):
    """ """
    assert type(dictionary) is dict, "The input must be a dictionary"

    metrics = dictionary.keys()
    sample_ids = dictionary.values()[0].keys()

    # gurantee that the data is ordered according to the sample identifiers
    _buffer = {}
    for m_key, m_value in dictionary.iteritems():
        for d_key, d_value in m_value.iteritems():
            if d_key not in _buffer:
                _buffer[d_key] = []
            _buffer[d_key].append(d_value)

    # unfold the dictionary and pack it in a list of lists
    data = []
    for sample_id in sample_ids:
        data.append(_buffer[sample_id])

    return metrics, sample_ids, array(data)
