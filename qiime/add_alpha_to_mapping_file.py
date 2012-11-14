#!/usr/bin/env python
# File created on 02 Nov 2012
from __future__ import division

__author__ = "Yoshiki Vazquez-Baeza"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Yoshiki Vazquez-Baeza", "Antonio Gonzalez-Pena"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Yoshiki Vazquez-Baeza"
__email__ = "yoshiki89@gmail.com"
__status__ = "Development"


from copy import deepcopy
from numpy import floor, ceil
from numpy import array, ndarray

def add_alpha_diversity_values_to_mapping_file(metrics, alpha_sample_ids,\
                                            alpha_data, mapping_file_headers,\
                                            mapping_file_data, bins,\
                                            missing_value_name='N/A'):
    """add 3 columns in the mapping file representing the alpha diversity data

    Inputs:
    metrics: list of alpha diversity metrics
    alpha_sample_ids: list of sample identifiers in the alpha diversity data
    alpha_data: alpha diversity data, as many columns as metrics and as many
    rows as elments in alpha_sample_ids
    mapping_file_headers: list of headers for the metadata mapping file
    mapping_file_data: metadata mapping file data
    bins: bins to classify the alpha diversity data
    missing_value_name: string to place for the sample ids in the mapping file
    but not in the alpha diversity data

    Output:
    mapping_file_headers: modified headers mapping file headers, including three
    new columns per metric i. e. for chao1 the new fields would be:
    'chao1_alpha', 'chao1_alpha_norm' and 'chao1_alpha_bin'
    mapping_file_data: extended data including the alpha diversity values,
    normalized values and bins

    """

    # data will be modified and returned so get your own copy
    new_mapping_file_data = deepcopy(mapping_file_data)
    new_mapping_file_headers = deepcopy(mapping_file_headers)

    for index, metric in enumerate(metrics):
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
            value.append(_get_level(norm_value, bins, 'bin'))

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
    """get a level from a value between 0 to 1; return a label or an integer

    Input:
    value: normalized value to assign a level to, must be between 0 and 1
    levels: number of in which a value can be assigned to
    prefix: a string prefix serves as a tagger

    Output:
    output: (str) if a prefix is provided a string of the form: prefix_2_of_5 is
    returned, where 2 is the level assigned to the passed value and 5 is the
    number of levels specified. (int) If a prefix is not provided the level of
    the value is returned as an integer

    """
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

def quantile(data, quantiles):
    """calculates quantiles of a dataset matching a given list of probabilities

    Input:
    data: 1-D list or numpy array with data to calculate the quantiles
    quantiles: list of probabilities, floating point values between 0 and 1

    Output:
    A list of elements drawn from 'data' that corresponding to the list of
    probabilities. This by default is using R. type 7 method for computation of
    the quantiles.
    """

    assert type(data) == list or type(data) == ndarray, "Data must be either"+\
        " a Python list or a NumPy 1-D array"
    assert type(quantiles) == list or type(quantiles) == ndarray, "Quantiles"+\
        "must be either a Python list or a NumPy 1-D array"
    assert all(map(lambda x: x>=0 and x<=1, quantiles)), "All the elements "+\
        "in the quantiles list must be greater than 0 and lower than one"

    # unless the user wanted, do not modify the data
    data = deepcopy(data)

    if type(data) != ndarray:
        data = array(data)

    data.sort()

    output = []

    # if needed different quantile methods could be used
    for one_quantile in quantiles:
        output.append(_quantile(data, one_quantile))

    return output

def _quantile(data, quantile):
    """gets a single quantile value for a dataset using R. type 7 method

    Input:
    data: sorted 1-d numpy array with float or int elements
    quantile: floating point value between 0 and 1

    Output:
    quantile value of data

    This function is based on cogent.maths.stats.util.NumbersI
    """
    index = quantile*(len(data)-1)
    bottom_index = int(floor(index))
    top_index = int(ceil(index))

    difference = index-bottom_index
    output = (1-difference)*data[bottom_index]+difference*data[top_index]

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
