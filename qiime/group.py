#!/usr/bin/env python

"""This module contains functions useful for obtaining groupings."""

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jai Ram Rideout",
               "Greg Caporaso",
               "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from collections import defaultdict
from functools import partial

import pandas as pd
import numpy as np

from qiime.stats import is_symmetric_and_hollow
from qiime.parse import group_by_field, parse_mapping_file
from qiime.filter import filter_mapping_file


def get_grouped_distances(dist_matrix_header, dist_matrix, mapping_header,
                          mapping, field, within=True,
                          suppress_symmetry_and_hollowness_check=False):
    """Returns a list of distance groupings for the specified field.

    The return value is a list that contains tuples of three elements: the
    first two elements are the field values being compared, and the third
    element is a list of the distances.

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.

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
        - suppress_symmetry_and_hollowness_check: By default, the input
          distance matrix will be checked for symmetry and hollowness. It is
          recommended to leave this check in place for safety, as the check
          is fairly fast. However, if you *know* you have a symmetric and
          hollow distance matrix, you can disable this check for small
          performance gains on extremely large distance matrices
    """
    _validate_input(dist_matrix_header, dist_matrix, mapping_header, mapping,
                    field)
    mapping_data = [mapping_header]
    mapping_data.extend(mapping)
    groups = group_by_field(mapping_data, field)
    return _get_groupings(dist_matrix_header, dist_matrix, groups, within,
                          suppress_symmetry_and_hollowness_check)


def get_all_grouped_distances(dist_matrix_header, dist_matrix, mapping_header,
                              mapping, field, within=True,
                              suppress_symmetry_and_hollowness_check=False):
    """Returns a list of distances for either samples within each of the
    field values or between each of the field values for the specified field.

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.

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
        - suppress_symmetry_and_hollowness_check: By default, the input
          distance matrix will be checked for symmetry and hollowness. It is
          recommended to leave this check in place for safety, as the check
          is fairly fast. However, if you *know* you have a symmetric and
          hollow distance matrix, you can disable this check for small
          performance gains on extremely large distance matrices
    """
    distances = get_grouped_distances(dist_matrix_header, dist_matrix,
                                      mapping_header, mapping, field, within,
                                      suppress_symmetry_and_hollowness_check)
    results = []
    for group in distances:
        for distance in group[2]:
            results.append(distance)
    return results


def get_field_state_comparisons(dist_matrix_header, dist_matrix,
                                mapping_header, mapping, field,
                                comparison_field_states,
                                suppress_symmetry_and_hollowness_check=False):
    """Returns a 2D dictionary relating distances between field states.

    The 2D dictionary is constructed such that each top-level key is a field
    state other than the field states in comparison_field_states. The
    second-level key is a field state from comparison_field_states, and the
    value at the (key, key) index is a list of distances between those two
    field states. Thus, given a field, this function will create comparisons
    between the specified comparison_field_states and all other field states.

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.

    Arguments:
        - dist_matrix_header: The distance matrix header, obtained from
                              parse.parse_distmat()
        - dist_matrix: The distance matrix, obtained from
                       parse.parse_distmat().
        - mapping_header: The mapping file header, obtained from
                          parse.parse_mapping_file()
        - mapping: The mapping file's contents, obtained from
                   parse.parse_mapping_file()
        - field: A field in the mapping file to do the comparisons on.
        - comparison_field_states: A list of strings specifying the field
          states to compare to all other field states. Cannot be an empty list.
        - suppress_symmetry_and_hollowness_check: By default, the input
          distance matrix will be checked for symmetry and hollowness. It is
          recommended to leave this check in place for safety, as the check
          is fairly fast. However, if you *know* you have a symmetric and
          hollow distance matrix, you can disable this check for small
          performance gains on extremely large distance matrices
    """
    _validate_input(dist_matrix_header, dist_matrix, mapping_header, mapping,
                    field)

    # avoid empty groups of distances
    mapping_header, mapping = filter_mapping_file(mapping, mapping_header,
                                                  dist_matrix_header)

    # Make sure each comparison group field state is in the specified field.
    if not comparison_field_states:
        raise ValueError("You must provide at least one field state to "
                         "compare to all of the other field states.")
    mapping_data = [mapping_header]
    mapping_data.extend(mapping)
    groups = group_by_field(mapping_data, field)
    for field_state in comparison_field_states:
        if field_state not in groups:
            raise ValueError("The comparison group field state '%s' is not in "
                             "the provided mapping file's field '%s'."
                             % (field_state, field))

    # Grab a list of all other field states (besides the ones in
    # comparison_field_states). These will be the field states that the states
    # in comparison_field_states will be compared against.
    field_states = [group for group in groups.keys()
                    if group not in comparison_field_states]

    # Get between distance groupings for the field of interest.
    between_groupings = get_grouped_distances(dist_matrix_header, dist_matrix,
                                              mapping_header, mapping, field, within=False,
                                              suppress_symmetry_and_hollowness_check=
                                              suppress_symmetry_and_hollowness_check)

    # Build up our 2D dictionary giving the distances between a field state and
    # a comparison group field state by filtering out the between_groupings
    # list to include only the comparisons that we want.
    result = {}
    for field_state in field_states:
        result[field_state] = {}
        for comp_field_state in comparison_field_states:
            result[field_state][comp_field_state] = []
            for group in between_groupings:
                if ((group[0] == field_state or group[1] == field_state)
                    and (group[0] == comp_field_state or
                         group[1] == comp_field_state)):
                    # We've found a group of distances between our comparison
                    # field state and the current field state, so keep the
                    # data.
                    result[field_state][comp_field_state] = group[2]
    return result


def get_ordered_coordinates(coordinate_header,
                            coordinate_matrix,
                            order,
                            strict=False):
    """ Return coordinate vectors in order

        coordinate_header: ids corresponding to vectors
         in coordinate_matrix (element 0 of output of
         qiime.parse.parse_coords)
        coordinate_matrix: the coordinate vectors (element 1 of
         output of qiime.parse.parse_coords)
        order: ordered ids from coordinate_header (usually sample
         ids) for coordinates that should be extracted
        strict: raise an error if an id from order is not present
         in coordinate_header (default: that id is ignored)

        The output of this function will be a tuple of the coordinate
         vectors corresponding to each id in order, and the id order:
         (ordered_coordinates, ordered_ids)
        Note that the output order can be a subset of the input order
         if some ids from order are not present in coordinate_header
         and strict == False.

        This function can be used in a way analogous to
         get_adjacent_distances to get a set of coordinates that
         might be connected by a line, for example.
    """
    ordered_coordinates = []
    ordered_ids = []
    for o in order:
        try:
            coordinate_idx = coordinate_header.index(o)
        except ValueError:
            if strict:
                raise ValueError(
                    "ID (%s) is not present in coordinate matrix" %
                    o)
            else:
                pass
        else:
            ordered_coordinates.append(coordinate_matrix[coordinate_idx])
            ordered_ids.append(o)
    return ordered_coordinates, ordered_ids


def get_adjacent_distances(dist_matrix_header,
                           dist_matrix,
                           sample_ids,
                           strict=False):
    """Return the distances between the adjacent sample_ids as a list

    dist_matrix_header: distance matrix headers, e.g. the output
        of qiime.parse.parse_distmat (element 0)
    dist_matrix: distance matrix, e.g., the output of
        qiime.parse.parse_distmat (element 1)
    sample_ids: a list of sample ids
    strict: boolean indicating whether to raise ValueError if a
        sample_id is not in dm (default: False; sample_ids not in
        dm are ignored)

    The output of this function will be a list of the distances
    between the adjacent sample_ids, and a list of the pair of sample ids
    corresponding to each distance. This could subsequently be used, for
    example, to plot unifrac distances between days in a timeseries, as
    d1 to d2, d2 to d3, d3 to d4, and so on. The list of pairs of sample
    ids are useful primarily in labeling axes when strict=False

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.

    """
    filtered_idx = []
    filtered_sids = []
    for sid in sample_ids:
        try:
            idx = dist_matrix_header.index(sid)
        except ValueError:
            if strict:
                raise ValueError(
                    "Sample ID (%s) is not present in distance matrix" %
                    sid)
            else:
                pass
        else:
            filtered_idx.append(idx)
            filtered_sids.append(sid)

    if len(filtered_idx) < 2:
        raise ValueError("At least two of your sample_ids must be present in the"
                         " distance matrix. %d are present." % len(filtered_idx))

    distance_results = []
    header_results = []
    for i in range(len(filtered_idx) - 1):
        distance_results.append(
            dist_matrix[filtered_idx[i]][filtered_idx[i + 1]])
        header_results.append(
            (filtered_sids[i], filtered_sids[i + 1]))
    return distance_results, header_results

def _group_by_sample_metadata(collapsed_md, sample_id_field="SampleID"):
    """Group sample identifiers by one or more metadata fields

    Parameters
    ----------
    collapsed_md : pd.DataFrame
        The result of collapsing a sample metadata DataFrame, for example with
        collapse_metadata.
    sample_id_field : str, optional
        The sample id field in the mapping_f.

    Returns
    -------
    dict
        Mapping of group id to set of input sample ids in that group.
    dict
        Mapping of input sample id to new group id.
    pd.DataFrame
        Sample metadata resulting from the collapse operation.

    Raises
    ------
    KeyError
        If sample_id_field or any of the collapse fields are not column headers
        in mapping_f.

    """
    new_index_to_group = {}
    old_index_to_new_index = {}
    for i in collapsed_md.index:
        old_indices = collapsed_md[sample_id_field][i]

        # this is a little ugly, but we need to handle single and multi-index
        # values here, and we always want to result to be a tuple
        if isinstance(i, tuple):
            new_index = i
        else:
            new_index = (i, )

        new_index_to_group[new_index] = set(old_indices)
        for old_index in old_indices:
            old_index_to_new_index[old_index] = new_index

    return new_index_to_group, old_index_to_new_index

def get_collapse_fns():
    """ Return lookup of functions that can be used with biom.Table.collapse
    """
    return {'median': _collapse_to_median,
            'first': _collapse_to_first,
            'random': _collapse_to_random,
            'sum': _collapse_to_sum,
            'mean': _collapse_to_mean}

def collapse_samples(table, mapping_f, collapse_fields, collapse_mode):
    """ Collapse samples in a biom table and sample metadata

    Parameters
    ----------
    table : biom.Table
        The biom table to be collapsed.
    mapping_f : file handle or filepath
        The sample metadata mapping file.
    collapse_fields : iterable
        The fields to combine when collapsing samples. For each sample in the
        mapping_f, the ordered values from these columns will be tuplized and
        used as the group identfier. Samples whose tuplized values in these
        fields are identical will be grouped.
    collapse_mode : str {sum, mean, median, random, first}
        The strategy to use for collapsing counts in the table.

    Returns
    -------
    biom.Table
        The collapsed biom table.
    pd.DataFrame
        Sample metadata resulting from the collapse operation.

    Raises
    ------
    KeyError
        If sample_id_field or any of the collapse fields are not column headers
        in mapping_f.

    """
    collapsed_metadata = _collapse_metadata(mapping_f,
                                            collapse_fields)

    new_index_to_group, old_index_to_new_index = \
        _group_by_sample_metadata(collapsed_metadata)
    partition_f = partial(_sample_id_from_group_id,
                          sid_to_group_id=old_index_to_new_index)

    collapse_fns = get_collapse_fns()
    try:
        collapse_f = collapse_fns[collapse_mode]
    except KeyError:
        raise KeyError(
         "Unknown collapse function %s. Valid choices are: "
         "%s." % (collapse_mode, ', '.join(collapse_fns.keys())))
    output_table = table.collapse(
        partition_f, collapse_f=collapse_f, norm=False, axis='sample')

    return collapsed_metadata, output_table

def mapping_lines_from_collapsed_df(collapsed_df):
    """ Formats a multi-index DataFrame as lines of a QIIME mapping file

    Parameters
    ----------
    collapsed_df : pd.DataFrame
        Sample metadata resulting from the collapse operation.

    Returns
    -------
    list of strings
        Lines representing the text of a QIIME mapping file.
    """
    lines = []
    lines.append('\t'.join(['#SampleID', 'original-sample-ids'] +\
                           list(collapsed_df.columns)[1:]))

    for r in collapsed_df.iterrows():
        # this is a little ugly, but we need to handle single and multi-index
        # values here
        if isinstance(r[0], tuple):
            new_idx = '.'.join(map(str, r[0]))
        else:
            new_idx = str(r[0])
        new_values = []
        for e in r[1]:
            if len(set(e)) == 1:
                # if all samples in the replicate group have the same
                # value for this column, just store that value
                new_values.append(str(e[0]))
            else:
                # if any samples in the replicate group differ in the value
                # in this column, store all of the values in the same order
                # as the ids in the new "original-sample-ids" column
                new_values.append('(%s)' % ', '.join(map(str,e)))
        lines.append('\t'.join([new_idx] + new_values))
    return lines

def _collapse_metadata(mapping_f, collapse_fields):
    """ Load a mapping file into a DataFrame and then collapse rows

    Parameters
    ----------
    mapping_f : file handle or filepath
        The sample metadata mapping file.
    collapse_fields : iterable
        The fields to combine when collapsing samples. For each sample in the
        mapping_f, the ordered values from these columns will be tuplized and
        used as the group identfier. Samples whose tuplized values in these
        fields are identical will be grouped.

    Returns
    -------
    pd.DataFrame
        Sample metadata resulting from the collapse operation.

    Raises
    ------
    KeyError
        If sample_id_field or any of the collapse fields are not column headers
        in mapping_f.

    """
    mapping_data, header, _ = parse_mapping_file(mapping_f)
    sample_md = pd.DataFrame(mapping_data, columns=header)
    grouped = sample_md.groupby(collapse_fields)
    collapsed_md = grouped.agg(lambda x: tuple(x))
    return collapsed_md

def _sample_id_from_group_id(id_, md, sid_to_group_id):
    try:
        group_id = sid_to_group_id[id_]
    except KeyError:
        raise KeyError("Sample id %s doesn't map to a group id." % id_)
    return '.'.join(map(str, group_id))

def _collapse_to_first(t, axis):
    return np.asarray([e[0] for e in t.iter_data(axis=axis, dense=True)])

def _collapse_to_median(t, axis):
    return np.asarray([np.median(e) for e in t.iter_data(axis=axis, dense=True)])

def _collapse_to_sum(t, axis):
    return np.asarray([np.sum(e) for e in t.iter_data(axis=axis)])

def _collapse_to_mean(t, axis):
    return np.asarray([np.mean(e) for e in t.iter_data(axis=axis)])

def _collapse_to_random(t, axis):
    if axis == 'sample':
        length = t.length("observation")
    elif axis == 'observation':
        length = t.length("sample")
    else:
        raise UnknownAxisError(axis)
    n = np.random.randint(length)
    return np.asarray([e[n] for e in t.iter_data(axis=axis, dense=True)])

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

    # check that we share sample identifiers between th mf and the dm
    if not set(zip(*mapping)[0]) & set(dist_matrix_header):
        raise ValueError('The mapping file does not share at least one sample'
                         ' with the distance matrix.')


def _get_indices(input_items, wanted_items):
    """Returns indices of the wanted items in the input items if present.

    input_items must be iterable, and wanted_items may be either a single value
    or a list. The return value will always be a list of indices, and an empty
    list if none were found. If wanted_items is a single string, it is treated
    as a scalar, not an iterable.
    """
    # Note: Some of this code is taken from Jeremy Widmann's
    # get_valid_indices() function, part of make_distance_histograms.py from QIIME 1.8.0.
    try:
        iter(input_items)
    except:
        raise ValueError("The input_items to search must be iterable.")
    try:
        len(wanted_items)
    except:
        # We have a scalar value, so put it in a list.
        wanted_items = [wanted_items]
    if isinstance(wanted_items, basestring):
        wanted_items = [wanted_items]

    return [input_items.index(item)
            for item in wanted_items if item in input_items]


def _get_groupings(dist_matrix_header, dist_matrix, groups, within=True,
                   suppress_symmetry_and_hollowness_check=False):
    """Returns a list of distance groupings.

    The return value is a list that contains tuples of three elements: the
    first two elements are the field values being compared, and the third
    element is a list of the distances.

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.

    Arguments:
        - dist_matrix_header: The distance matrix header.
        - dist_matrix: The distance matrix.
        - groups: A dictionary mapping field value to sample IDs, obtained by
                  calling group_by_field().
        - within: If True, distances are grouped within a field value. If
          False, distances are grouped between field values.
        - suppress_symmetry_and_hollowness_check: By default, the input
          distance matrix will be checked for symmetry and hollowness. It is
          recommended to leave this check in place for safety, as the check
          is fairly fast. However, if you *know* you have a symmetric and
          hollow distance matrix, you can disable this check for small
          performance gains on extremely large distance matrices

    If within is True, the zeros along the diagonal of the distance matrix are
    omitted.
    """
    # Note: Much of this code is taken from Jeremy Widmann's
    # distances_by_groups() function, part of make_distance_histograms.py from QIIME 1.8.0.
    if not suppress_symmetry_and_hollowness_check:
        if not is_symmetric_and_hollow(dist_matrix):
            raise ValueError("The distance matrix must be symmetric and "
                             "hollow.")
    result = []
    group_items = groups.items()

    for i, (row_group, row_ids) in enumerate(group_items):
        row_indices = _get_indices(dist_matrix_header, row_ids)
        if within:
            # Handle the case where indices are the same so we need to omit
            # the diagonal.
            block = dist_matrix[row_indices][:, row_indices]

            size = len(row_indices)
            indices = []
            for i in range(size):
                for j in range(i, size):
                    if i != j:
                        indices.append(block[i][j])
            if indices:
                result.append((row_group, row_group, indices))
        else:
            # Handle the case where indices are separate: just return blocks.
            for j in range(i + 1, len(groups)):
                col_group, col_ids = group_items[j]
                col_indices = _get_indices(dist_matrix_header, col_ids)
                vals = dist_matrix[row_indices][:, col_indices]

                # Flatten the array into a single-level list.
                vals = map(None, vals.flat)
                if vals:
                    result.append((row_group, col_group, vals))
    return result


def extract_per_individual_states_from_sample_metadata(
        sample_metadata,
        state_category,
        state_values,
        individual_identifier_category,
        filter_missing_data=True):
    """
    sample_metadata : 2d dictionary mapping sample ids to metadata (as
     returned from qiime.parse.parse_mapping_file_to_dict)
    state_category: metadata category name describing state of interest
     (usually something like 'TreatmentState') as a string
    state_values: ordered list of values of interest in the state_category
     metadata entry (usually something like ['PreTreatment','PostTreatment'])
    individual_identifier_category: metadata category name describing the
     individual (usually something like 'PersonalID') as a string
    filter_missing_data: if True, an individual is excluded
     from the result object if any of it's values are None. This can occur
     when there is no sample for one or more of the state values for an
     individual. This is True by default.

    returns {'individual-identifier':
               [sample-id-at-state-value1,
                sample-id-at-state-value2,
                sample-id-at-state-value3,
                ...],
              ...
             }
    """
    # prep the result object, which will be a dict of lists
    len_state_values = len(state_values)

    def inner_dict_constructor():
        return [None] * len_state_values
    results = defaultdict(inner_dict_constructor)

    for sample_id, metadata in sample_metadata.items():
        try:
            individual_id = metadata[individual_identifier_category]
        except KeyError:
            raise KeyError("%s is not a sample metadata category." %
                           individual_identifier_category)
        try:
            state_value = metadata[state_category]
        except KeyError:
            raise KeyError("%s is not a sample metadata category." %
                           state_category)

        try:
            state_index = state_values.index(state_value)
        except ValueError:
            # hit a state that is in the mapping file but not in
            # state_values - this is silently ignored
            continue

        results[individual_id][state_index] = sample_id

    if filter_missing_data:
        # delete individual results if sample ids corresponding to
        # any of the states are missing
        for individual_id, sample_ids in results.items():
            if None in sample_ids:
                del results[individual_id]
    return results


def extract_per_individual_state_metadatum_from_sample_metadata(
        sample_metadata,
        state_category,
        state_values,
        individual_identifier_category,
        metadata_category,
        process_f=float):
    """
    sample_metadata : 2d dictionary mapping sample ids to metadata (as
     returned from qiime.parse.parse_mapping_file_to_dict)
    state_category: metadata category name describing state of interest
     (usually something like 'TreatmentState') as a string
    state_values: ordered list of values of interest in the state_category
     metadata entry (usually something like ['PreTreatment','PostTreatment'])
    individual_identifier_category: metadata category name describing the
     individual (usually something like 'PersonalID') as a string
    metadata_category: metadata category to extract from sample_metadata
    process_f: function to apply to metadata values (default: float)

    returns {'individual-identifier':
               [state-1-metadata-value,
                state-2-metadata-value,
                ...],
              ...
             }
    """
    per_individual_states = extract_per_individual_states_from_sample_metadata(
        sample_metadata,
        state_category,
        state_values,
        individual_identifier_category,
        filter_missing_data=True)

    results = {}
    for individual_id, sample_ids in per_individual_states.items():
        per_state_metadata_values = []
        for sample_id in sample_ids:
            try:
                sample_metadata_value = sample_metadata[
                    sample_id][
                    metadata_category]
            except KeyError:
                raise KeyError(
                    "%s is not a sample metadata category." %
                    metadata_category)
            try:
                v = process_f(sample_metadata_value)
            except ValueError as e:
                v = None
            per_state_metadata_values.append(v)
        results[individual_id] = per_state_metadata_values
    return results


def extract_per_individual_state_metadata_from_sample_metadata(
        sample_metadata,
        state_category,
        state_values,
        individual_identifier_category,
        metadata_categories,
        process_f=float):
    """
    sample_metadata : 2d dictionary mapping sample ids to metadata (as
     returned from qiime.parse.parse_mapping_file_to_dict)
    state_category: metadata category name describing state of interest
     (usually something like 'TreatmentState') as a string
    state_values: ordered list of values of interest in the state_category
     metadata entry (usually something like ['PreTreatment','PostTreatment'])
    individual_identifier_category: metadata category name describing the
     individual (usually something like 'PersonalID') as a string
    metadata_categories: metadata categories to extract from sample_metadata
    process_f: function to apply to metadata values (default: float)

    returns {'metadata-category-1':
              {'individual-identifier-1':
               [difference-in-metadata-value-bw-states-2-and-1,
                difference-in-metadata-value-bw-states-3-and-2,
                ...],
               'individual-identifier-2:
               [difference-in-metadata-value-bw-states-2-and-1,
                difference-in-metadata-value-bw-states-3-and-2,
                ...],
               }
              ...
              }
    """
    results = {}
    for metadata_category in metadata_categories:
        results[metadata_category] = \
            extract_per_individual_state_metadatum_from_sample_metadata(
                sample_metadata,
                state_category,
                state_values,
                individual_identifier_category,
                metadata_category,
                process_f)
    return results


def extract_per_individual_state_metadata_from_sample_metadata_and_biom(
        sample_metadata,
        biom_table,
        state_category,
        state_values,
        individual_identifier_category,
        observation_ids=None):
    """
    sample_metadata : 2d dictionary mapping sample ids to metadata (as
     returned from qiime.parse.parse_mapping_file_to_dict)
    biom_table: biom table object containing observation counts for
     samples in sample_metadata
    state_category: metadata category name describing state of interest
     (usually something like 'TreatmentState') as a string
    state_values: ordered list of values of interest in the state_category
     metadata entry (usually something like ['PreTreatment','PostTreatment'])
    individual_identifier_category: metadata category name describing the
     individual (usually something like 'PersonalID') as a string
    observation_ids: observations (usually OTUs) to extract from biom_table
     (default is all)

    returns {'otu1':
              {'individual-identifier-1:
               [difference-in-otu1-abundance-bw-states-2-and-1,
                difference-in-otu1-abundance-bw-states-3-and-2,
                ...],
               'individual-identifier-2:
               [difference-in-otu1-abundance-bw-states-2-and-1,
                difference-in-otu1-abundance-bw-states-3-and-2,
                ...],
               }
              ...
              }
    """
    per_individual_states = extract_per_individual_states_from_sample_metadata(
        sample_metadata,
        state_category,
        state_values,
        individual_identifier_category,
        filter_missing_data=True)
    results = {}
    if observation_ids is None:
        observation_ids = biom_table.ids(axis='observation')
    for observation_id in observation_ids:
        observation_data = biom_table.data(observation_id, 'observation')
        results[observation_id] = {}
        for individual_id, sample_ids in per_individual_states.items():
            per_state_metadata_values = []
            for sample_id in sample_ids:
                sample_index = biom_table.index(sample_id, 'sample')
                per_state_metadata_values.append(
                    observation_data[sample_index])
            results[observation_id][individual_id] = per_state_metadata_values
    return results
