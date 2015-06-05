#!/usr/bin/env python
# File created on 02 Nov 2012
from __future__ import division

__author__ = "Yoshiki Vazquez-Baeza"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Yoshiki Vazquez-Baeza", "Antonio Gonzalez-Pena"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Yoshiki Vazquez-Baeza"
__email__ = "yoshiki89@gmail.com"


from copy import deepcopy
from qiime.stats import quantile
from numpy import searchsorted, array
from qiime.parse import parse_rarefaction


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
    norm = lambda x, x_min, x_max: (x - x_min) / (x_max - x_min)

    # data will be modified and returned so get your own copy
    new_mapping_file_data = deepcopy(mapping_file_data)
    new_mapping_file_headers = deepcopy(mapping_file_headers)

    # regular levels assigned based on equally spaced bins
    overall_probs = [i / bins for i in range(1, bins)]

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
        # metric being used; hence the calculation and normalization of the
        # data
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

    if value > 1 or value < 0:
        raise ValueError("Encountered invalid normalized alpha diversity value %s. "
                         "Normalized values must be between 0 and 1." % value)

    check = [i for i in range(0, len(levels)) if levels[i] == value]

    # apply a special rule for the values that are equal to an edge
    if len(check):
        value_level = check[0] + 2
    # if it is not a special case just use searchsorted
    else:
        value_level = searchsorted(levels, value) + 1

    if prefix is not None:
        output = '{0}_{1}_of_{2}'.format(prefix, value_level, len(levels) + 1)
    else:
        output = value_level

    return output


def mean_alpha(alpha_dict, depth):
    """mean collated alpha diversity data at a given depth

    Input:
    alpha_dict: dictionary where the values are the lines of a collated alpha
    diversity data files and the keys are the names of each of these files with
    no extension, this name is usually the metric used to compute the alpha
    diversity.
    depth: selected depth to mean the computed alpha diversity values for the
    alpha_dict data. If None is passed, the highest depth will be used.

    Output:
    metrics: list of metric names i. e. the name of each collated alpha div file
    sample_ids: list of sample identifiers represented
    data: a list of lists with the mean of alpha diversity data at a given
    depth for the different metrics, each column is a different metric.
    """

    assert isinstance(alpha_dict, dict), "Input data must be a dictionary"
    assert depth is None or (depth >= 0 and isinstance(depth, int)), "The " +\
        "specified depth must be a positive integer."

    metrics = []
    sample_ids = []
    data = []

    for key, value in alpha_dict.iteritems():
        identifiers, _, _, rarefaction_data = parse_rarefaction(value)

        # if depth is specified as None use the highest available, retrieve it
        # on a per file basis so you make sure the value exists for all files
        if depth is None:
            _depth = int(max([row[0] for row in rarefaction_data]))
        else:
            _depth = depth
        metrics.append('{0}_even_{1}'.format(key, _depth))

        # check there are elements with the desired rarefaction depth
        if sum([1 for row in rarefaction_data if row[0] == _depth]) == 0:
            # get a sorted list of strings with the available rarefaction
            # depths
            available_rarefaction_depths = map(str, sorted(list(set([row[0] for
                                                                     row in rarefaction_data]))))
            raise ValueError("The depth %d does not exist in the collated "
                             "alpha diversity file for the metric: %s. The available depths "
                             "are: %s." % (_depth, key, ', '.join(available_rarefaction_depths)))

        # check all the files have the same sample ids in the same order
        if sample_ids:
            if not sample_ids == identifiers[3:]:
                raise ValueError("Non-matching sample ids were found in the "
                                 "collated alpha diversity files. Make sure all the files "
                                 "contain data for the same samples.")
        else:
            sample_ids = identifiers[3:]

        # find all the data at the desired depth and get the mean values, remove
        # the first two elements ([depth, iteration]) as those are not needed
        data.append(array([row[2:] for row in rarefaction_data if
                           row[0] == _depth]).mean(axis=0))

    # transpose the data to match the formatting of non-collated alpha div data
    data = array(data).T.tolist()

    return metrics, sample_ids, data
