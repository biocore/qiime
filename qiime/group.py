#!/usr/bin/env python

"""This module contains functions useful for obtaining groupings."""

__author__ = "Jai Rideout"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jai Rideout", "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Jai Rideout"
__email__ = "jr378@nau.edu"
__status__ = "Development"

from numpy import array
from qiime.parse import group_by_field

def get_grouped_distances(dist_matrix_header, dist_matrix, mapping_header,
                          mapping, field, within=True):
    """Returns a list of distance groupings for the specified field.

    The return value is a list that contains tuples of three elements: the
    first two elements are the field values being compared, and the third
    element is a list of the distances.

    Arguments:
        - dist_matrix_header: The distance matrix header, obtained from
                              parse.parse_distmat()
        - dist_matrix: The distance matrix, obtained from
                       parse.parse_distmat().
        - mapping_header: The mapping file header, obtained from
                          parse.parse_mapping_file()
        - mapping: The mapping file's contents, obtained from
                   parse.parse_mapping_file()
        - field: A field in the mapping file to do the grouping on.
        - within: If True, distances are grouped within a field value. If
          False, distances are grouped between field values.
    """
    _validate_input(dist_matrix_header, dist_matrix, mapping_header, mapping,
                    field)
    mapping_data = [mapping_header]
    mapping_data.extend(mapping)
    groups = group_by_field(mapping_data, field)
    return _get_groupings(dist_matrix_header, dist_matrix, groups, within)

def get_all_grouped_distances(dist_matrix_header, dist_matrix, mapping_header,
                              mapping, field, within=True):
    """Returns a list of distances for either samples within each of the
    field values or between each of the field values for the specified field.

    Arguments:
        - dist_matrix_header: The distance matrix header, obtained from
                              parse.parse_distmat()
        - dist_matrix: The distance matrix, obtained from
                       parse.parse_distmat().
        - mapping_header: The mapping file header, obtained from
                          parse.parse_mapping_file()
        - mapping: The mapping file's contents, obtained from
                   parse.parse_mapping_file()
        - field: A field in the mapping file to do the grouping on.
        - within: If True, distances are grouped within a field value. If
          False, distances are grouped between field values.
    """
    distances = get_grouped_distances(dist_matrix_header, dist_matrix,
                                      mapping_header, mapping, field, within)
    results = []
    for group in distances:
        for distance in group[2]:
            results.append(distance)
    return results

def _validate_input(dist_matrix_header, dist_matrix, mapping_header, mapping,
                    field):
    """Validates the input data to make sure it can be used and makes sense.

    The headers, distance matrix, and mapping input should be iterable, and all
    data should not be None. The field must exist in the mapping header.
    """
    if (dist_matrix_header is None or dist_matrix is None or mapping_header is
        None or mapping is None or field is None):
        raise ValueError("The input(s) cannot be 'None'.")

    # Make sure the appropriate input is iterable.
    for input_arg in (dist_matrix_header, dist_matrix, mapping_header,
                      mapping):
        try:
            iter(input_arg)
        except:
            raise ValueError("The headers, distance matrix, and mapping data "
                             "must be iterable.")

    # The field must be a string.
    if not isinstance(field, str):
        raise ValueError("The field must be a string.")

    # Make sure the field is in the mapping header.
    if field not in mapping_header:
        raise ValueError("The field '%s' is not in the mapping file header."
                         % field)

def _get_indices(input_items, wanted_items):
    """Returns indices of the wanted items in the input items if present."""
    # Note: This code is taken from Jeremy Widmann's get_valid_indices()
    # function, part of make_distance_histograms.py.
    return [input_items.index(item)
            for item in wanted_items if item in input_items]

def _get_groupings(dist_matrix_header, dist_matrix, groups, within=True):
    """Returns a list of distance groupings.

    The return value is a list that contains tuples of three elements: the
    first two elements are the field values being compared, and the third
    element is a list of the distances.

    Arguments:
        - dist_matrix_header: The distance matrix header.
        - dist_matrix: The distance matrix.
        - groups: A dictionary mapping field value to sample IDs, obtained by
                  calling group_by_field().
        - within: If True, distances are grouped within a field value. If
          False, distances are grouped between field values.
    
    If within is True, the zeros along the diagonal of the distance matrix are
    omitted.
    """
    # Note: Much of this code is taken from Jeremy Widmann's
    # distances_by_groups() function, part of make_distance_histograms.py.
    result = []
    group_items = groups.items()

    for i, (row_group, row_ids) in enumerate(group_items):
        row_indices = _get_indices(dist_matrix_header, row_ids)
        if within:
            # Handle the case where indices are the same so we need to omit
            # the diagonal.
            block = dist_matrix[row_indices][:,row_indices]

            size = len(row_indices)
            indices = []
            for i in range(size):
                for j in range(i,size):
                    if i != j:
                        indices.append(block[i][j])
            if indices:
                result.append((row_group, row_group, indices))
        else:
            # Handle the case where indices are separate: just return blocks.
            for j in range(i+1, len(groups)):
                col_group, col_ids = group_items[j]
                col_indices = _get_indices(dist_matrix_header, col_ids)
                vals = dist_matrix[row_indices][:,col_indices]

                # Flatten the array into a single-level list.
                vals = map(None, vals.flat)
                if vals:
                    result.append((row_group, col_group, vals))
    return result
