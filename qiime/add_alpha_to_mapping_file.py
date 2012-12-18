#!/usr/bin/env python
# File created on 02 Nov 2012
from __future__ import division

__author__ = "Yoshiki Vazquez-Baeza"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Yoshiki Vazquez-Baeza", "Antonio Gonzalez-Pena"]
__license__ = "GPL"
__version__ = "1.6.0"
__maintainer__ = "Yoshiki Vazquez-Baeza"
__email__ = "yoshiki89@gmail.com"
__status__ = "Release"


from copy import deepcopy
from numpy import searchsorted
from qiime.stats import quantile

def add_alpha_diversity_values_to_mapping_file(metrics, alpha_sample_ids,
                                            alpha_data, mapping_file_headers,
                                            mapping_file_data, bins, method='equal',
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
    method: binning method selection, the options are 'equal' and 'quantile'.
    'equal', will get you equally spaced limits and 'quantile' will assign the
    limits using quantiles, using the selected number of bins.

    missing_value_name: string to place for the sample ids in the mapping file
    but not in the alpha diversity data

    Output:
    mapping_file_headers: modified headers mapping file headers, including three
    new columns per metric i. e. for chao1 the new fields would be:
    'chao1_alpha', 'chao1_alpha_norm' and 'chao1_alpha_bin'
    mapping_file_data: extended data including the alpha diversity values,
    normalized values and bins

    """
    norm = lambda x, x_min, x_max: (x-x_min)/(x_max-x_min)

    # data will be modified and returned so get your own copy
    new_mapping_file_data = deepcopy(mapping_file_data)
    new_mapping_file_headers = deepcopy(mapping_file_headers)

    # regular levels assigned based on equally spaced bins
    overall_probs = [i/bins for i in range(1, bins)]

    # if we are using the average method the levels are equal to the probs list
    if method == 'equal':
        levels = overall_probs

    for index, metric in enumerate(metrics):
        # get the alpha diversity value for the metric being evaluated
        data = [[row[index]] for row in alpha_data]
        metric_max = max(data)[0]
        metric_min = min(data)[0]

        # add headers for each metric
        new_mapping_file_headers.append('{0}_alpha'.format(metric))
        new_mapping_file_headers.append('{0}_normalized_alpha'.format(metric))
        new_mapping_file_headers.append('{0}_alpha_label'.format(metric))

        # when using the quantile method the levels change depending on the
        # metric being used; hence the calculation and normalization of the data
        if method == 'quantile':
            levels = quantile([norm(element[0], metric_min, metric_max)
                for element in data], overall_probs)

        # get the normalized value of diversity and the tag for each value
        for value in data:
            norm_value = norm(value[0], metric_min, metric_max)
            value.append(norm_value)
            value.append(_get_level(norm_value, levels, 'bin'))

        # iterate using the mapping file instead of using the alpha diversity
        # data because more often that you will like you will have more samples
        # in the mapping file than in your downstream analysis data
        for row in new_mapping_file_data:
            try:
                data_index = alpha_sample_ids.index(row[0])

                # data fields should be strings
                row.extend(map(str, data[data_index]))
            except ValueError:
                row.extend([missing_value_name, missing_value_name,
                    missing_value_name])

    return new_mapping_file_data, new_mapping_file_headers

def _get_level(value, levels, prefix=None):
    """accommodate a value into the 'levels' list; return a string or an integer

    Input:
    value: normalized value to assign a level to, must be between 0 and 1
    levels: sorted edges to check at which level does 'value' fall into place

    Output:
    output: (str) if a prefix is provided a string of the form: prefix_2_of_5 is
    returned, where 2 is the level assigned to the passed value and 5 is the
    number of levels specified. (int) If a prefix is not provided the level of
    the value is returned as an integer

    """
    assert value <= 1 and value >= 0, "The value must be between 0 and 1"

    check = [i for i in range(0, len(levels)) if levels[i] == value]

    # apply a special rule for the values that are equal to an edge
    if len(check):
        value_level = check[0] + 2
    # if it is not a special case just use searchsorted
    else:
        value_level = searchsorted(levels, value)+1

    if prefix != None:
        output = '{0}_{1}_of_{2}'.format(prefix, value_level, len(levels)+1)
    else:
        output = value_level

    return output
