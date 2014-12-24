#!/usr/bin/env python

"""Tests public and private functions in the group module."""

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jai Ram Rideout",
               "Greg Caporaso",
               "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.9.0-rc1"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from unittest import TestCase, main
from StringIO import StringIO

import numpy as np
from biom import Table
from biom.parse import parse_biom_table
from biom.exception import UnknownIDError
from numpy import array, matrix
from numpy.testing import assert_almost_equal
from StringIO import StringIO

from qiime.parse import (parse_mapping_file, parse_distmat,
                         group_by_field, parse_coords,
                         parse_mapping_file_to_dict)
from qiime.group import (get_grouped_distances, get_all_grouped_distances,
                         get_field_state_comparisons, _get_indices,
                         _get_groupings, _validate_input,
                         get_adjacent_distances, get_ordered_coordinates,
                         extract_per_individual_states_from_sample_metadata,
                         extract_per_individual_state_metadatum_from_sample_metadata,
                         extract_per_individual_state_metadata_from_sample_metadata_and_biom,
                         _group_by_sample_metadata, _sample_id_from_group_id,
                         collapse_samples, mapping_lines_from_collapsed_df,
                         _collapse_to_first, _collapse_to_median,
                         _collapse_to_random, _collapse_metadata,
                         _collapse_to_sum, _collapse_to_mean, get_collapse_fns)


class GroupTests(TestCase):

    """Tests of the group module."""

    def setUp(self):
        """Create some data to be used in the tests."""
        # Create the mapping file/distance matrix combo from the overview
        # tutorial.
        self.dist_matrix_string = ["\tPC.354\tPC.355\tPC.356\tPC.481\tPC.593\
                                    \tPC.607\tPC.634\tPC.635\tPC.636",
                                   "PC.354\t0.0\t0.625\t0.623\t0.61\t0.577\
                                    \t0.729\t0.8\t0.721\t0.765",
                                   "PC.355\t0.625\t0.0\t0.615\t0.642\t0.673\
                                    \t0.776\t0.744\t0.749\t0.677",
                                   "PC.356\t0.623\t0.615\t0.0\t0.682\t0.737\
                                    \t0.734\t0.777\t0.733\t0.724",
                                   "PC.481\t0.61\t0.642\t0.682\t0.0\t0.704\
                                    \t0.696\t0.675\t0.654\t0.696",
                                   "PC.593\t0.577\t0.673\t0.737\t0.704\t0.0\
                                    \t0.731\t0.758\t0.738\t0.737",
                                   "PC.607\t0.729\t0.776\t0.734\t0.696\t0.731\
                                    \t0.0\t0.718\t0.666\t0.727",
                                   "PC.634\t0.8\t0.744\t0.777\t0.675\t0.758\
                                    \t0.718\t0.0\t0.6\t0.578",
                                   "PC.635\t0.721\t0.749\t0.733\t0.654\t0.738\
                                    \t0.666\t0.6\t0.0\t0.623",
                                   "PC.636\t0.765\t0.677\t0.724\t0.696\t0.737\
                                    \t0.727\t0.578\t0.623\t0.0"]

        self.mapping_string = ["#SampleID\tBarcodeSequence\tTreatment\tDOB",
                               "PC.354\tAGCACGAGCCTA\tControl\t20061218",
                               "PC.355\tAACTCGTCGATG\tControl\t20061218",
                               "PC.356\tACAGACCACTCA\tControl\t20061126",
                               "PC.481\tACCAGCGACTAG\tControl\t20070314",
                               "PC.593\tAGCAGCACTTGT\tControl\t20071210",
                               "PC.607\tAACTGTGCGTAC\tFast\t20071112",
                               "PC.634\tACAGAGTCGGCT\tFast\t20080116",
                               "PC.635\tACCGCAGAGTCA\tFast\t20080116",
                               "PC.636\tACGGTGAGTGTC\tFast\t20080116"]

        # Field to test on. Field values are either "Control" or "Fast".
        self.field = 'Treatment'

        # Create a tiny distancy matrix/mapping file with a single sample for
        # additional testing.
        self.tiny_dist_matrix_string = ["\tSamp.1", "Samp.1\t0"]
        self.tiny_mapping_string = ["#SampleID\tBarcodeSequence\tSampleField",
                                    "Samp.1\tAGCACGAGCCTA\tSampleFieldState1"]
        self.tiny_field = 'SampleField'

        self.small_dist_matrix_string = ["\tSamp.1\tSamp.2", "Samp.1\t0\t0.5",
                                         "Samp.2\t0.5\t0"]
        self.small_mapping_string = ["#SampleID\tBarcodeSequence\tSampleField",
                                     "Samp.1\tAGCACGAGCCTA\tSampleFieldState1",
                                     "Samp.2\tAGCACGAGCCTG\tSampleFieldState2"]
        self.small_field = 'SampleField'

        # Parse mapping "files" (faked here).
        self.mapping, self.mapping_header, self.comments = parse_mapping_file(
            self.mapping_string)
        mapping_data = [self.mapping_header]
        mapping_data.extend(self.mapping)
        self.groups = group_by_field(mapping_data, self.field)

        self.tiny_mapping, self.tiny_mapping_header, self.tiny_comments = \
            parse_mapping_file(self.tiny_mapping_string)
        tiny_mapping_data = [self.tiny_mapping_header]
        tiny_mapping_data.extend(self.tiny_mapping)
        self.tiny_groups = group_by_field(tiny_mapping_data, self.tiny_field)

        self.small_mapping, self.small_mapping_header, self.small_comments = \
            parse_mapping_file(self.small_mapping_string)
        small_mapping_data = [self.small_mapping_header]
        small_mapping_data.extend(self.small_mapping)
        self.small_groups = group_by_field(small_mapping_data,
                                           self.small_field)

        # Parse distance matrix "files" (faked here).
        self.dist_matrix_header, self.dist_matrix = parse_distmat(
            self.dist_matrix_string)

        self.tiny_dist_matrix_header, self.tiny_dist_matrix = parse_distmat(
            self.tiny_dist_matrix_string)

        self.small_dist_matrix_header, self.small_dist_matrix = parse_distmat(
            self.small_dist_matrix_string)

        # extract_per_individual* input data
        self.individual_states_and_responses_map_f1 = \
            parse_mapping_file_to_dict(
                individual_states_and_responses_map_f1.split('\n'))[0]
        self.individual_states_and_responses_map_f2 = \
            parse_mapping_file_to_dict(
                individual_states_and_responses_map_f2.split('\n'))[0]
        self.paired_difference_biom1 = \
            parse_biom_table(paired_difference_biom_f1.split('\n'))

        self._group_by_sample_metadata_map_f1 = _group_by_sample_metadata_map_f1

    def test_get_grouped_distances_within(self):
        """get_grouped_distances() should return a list of within distance
        groupings."""
        groupings = get_grouped_distances(self.dist_matrix_header,
                                          self.dist_matrix, self.mapping_header, self.mapping,
                                          self.field, within=True)
        expected = [
            ('Control', 'Control', [0.625, 0.623, 0.60999999999999999,
                                    0.57699999999999996, 0.61499999999999999,
                                    0.64200000000000002, 0.67300000000000004,
                                    0.68200000000000005, 0.73699999999999999,
                                    0.70399999999999996]),
            ('Fast', 'Fast', [0.71799999999999997, 0.66600000000000004,
                              0.72699999999999998, 0.59999999999999998,
                              0.57799999999999996, 0.623])]
        self.assertEqual(groupings, expected)

    def test_get_grouped_distances_between(self):
        """get_grouped_distances() should return a list of between distance
        groupings."""
        groupings = get_grouped_distances(self.dist_matrix_header,
                                          self.dist_matrix, self.mapping_header, self.mapping,
                                          self.field, within=False)
        expected = [
            ('Control', 'Fast', [0.72899999999999998, 0.80000000000000004,
                                 0.72099999999999997, 0.76500000000000001,
                                 0.77600000000000002, 0.74399999999999999,
                                 0.749, 0.67700000000000005,
                                 0.73399999999999999, 0.77700000000000002,
                                 0.73299999999999998, 0.72399999999999998,
                                 0.69599999999999995, 0.67500000000000004,
                                 0.65400000000000003, 0.69599999999999995,
                                 0.73099999999999998, 0.75800000000000001,
                                 0.73799999999999999, 0.73699999999999999])]
        self.assertEqual(groupings, expected)

    def test_get_all_grouped_distances_within(self):
        """get_all_grouped_distances() should return a list of distances for
        all samples with the same field value."""
        groupings = get_all_grouped_distances(self.dist_matrix_header,
                                              self.dist_matrix, self.mapping_header, self.mapping,
                                              self.field, within=True)
        expected = [0.625, 0.623, 0.60999999999999999, 0.57699999999999996,
                    0.61499999999999999, 0.64200000000000002,
                    0.67300000000000004, 0.68200000000000005,
                    0.73699999999999999, 0.70399999999999996,
                    0.71799999999999997, 0.66600000000000004,
                    0.72699999999999998, 0.59999999999999998,
                    0.57799999999999996, 0.623]
        self.assertEqual(groupings, expected)

    def test_get_all_grouped_distances_between(self):
        """get_all_grouped_distances() should return a list of distances
        between samples of all different field values."""
        groupings = get_all_grouped_distances(self.dist_matrix_header,
                                              self.dist_matrix, self.mapping_header, self.mapping,
                                              self.field, within=False)
        expected = [0.72899999999999998, 0.80000000000000004,
                    0.72099999999999997, 0.76500000000000001,
                    0.77600000000000002, 0.74399999999999999, 0.749,
                    0.67700000000000005, 0.73399999999999999,
                    0.77700000000000002, 0.73299999999999998,
                    0.72399999999999998, 0.69599999999999995,
                    0.67500000000000004, 0.65400000000000003,
                    0.69599999999999995, 0.73099999999999998,
                    0.75800000000000001, 0.73799999999999999,
                    0.73699999999999999]
        self.assertEqual(groupings, expected)

    def test_get_field_state_comparisons(self):
        """get_field_state_comparisons() should return a 2D dictionary of
        distances between a field state and its comparison field states."""
        comparison_groupings = get_field_state_comparisons(
            self.dist_matrix_header, self.dist_matrix, self.mapping_header,
            self.mapping, self.field, ['Control'])
        expected = {'Fast': {'Control': [0.72899999999999998,
                                         0.80000000000000004, 0.72099999999999997, 0.76500000000000001,
                                         0.77600000000000002, 0.74399999999999999, 0.749,
                                         0.67700000000000005, 0.73399999999999999, 0.77700000000000002,
                                         0.73299999999999998, 0.72399999999999998, 0.69599999999999995,
                                         0.67500000000000004, 0.65400000000000003, 0.69599999999999995,
                                         0.73099999999999998, 0.75800000000000001, 0.73799999999999999,
                                         0.73699999999999999]}}
        self.assertDictEqual(comparison_groupings, expected)

        comparison_groupings = get_field_state_comparisons(
            self.dist_matrix_header, self.dist_matrix, self.mapping_header,
            self.mapping, self.field, ['Fast'])
        expected = {'Control': {'Fast': [0.72899999999999998,
                                         0.80000000000000004, 0.72099999999999997, 0.76500000000000001,
                                         0.77600000000000002, 0.74399999999999999, 0.749,
                                         0.67700000000000005, 0.73399999999999999, 0.77700000000000002,
                                         0.73299999999999998, 0.72399999999999998, 0.69599999999999995,
                                         0.67500000000000004, 0.65400000000000003, 0.69599999999999995,
                                         0.73099999999999998, 0.75800000000000001, 0.73799999999999999,
                                         0.73699999999999999]}}
        self.assertDictEqual(comparison_groupings, expected)

    def test_get_field_state_comparisons_extra_samples(self):
        self.small_mapping.append(['a', 'ACTCGAGGACT', 'xx'])
        comparison_groupings = get_field_state_comparisons(
            self.small_dist_matrix_header, self.small_dist_matrix,
            self.small_mapping_header, self.small_mapping,
            self.small_field, ['SampleFieldState1'])
        expected = {'SampleFieldState2': {'SampleFieldState1': [0.5]}}
        self.assertDictEqual(comparison_groupings, expected)

    def test_get_field_state_comparisons_small(self):
        """get_field_state_comparisons() should return a 2D dictionary of
        distances between a field state and its comparison field states."""
        comparison_groupings = get_field_state_comparisons(
            self.small_dist_matrix_header, self.small_dist_matrix,
            self.small_mapping_header, self.small_mapping,
            self.small_field, ['SampleFieldState1'])
        expected = {'SampleFieldState2': {'SampleFieldState1': [0.5]}}
        self.assertDictEqual(comparison_groupings, expected)

    def test_get_field_state_comparisons_tiny(self):
        """get_field_state_comparisons() should return an empty dictionary."""
        comparison_groupings = get_field_state_comparisons(
            self.tiny_dist_matrix_header, self.tiny_dist_matrix,
            self.tiny_mapping_header, self.tiny_mapping, self.tiny_field,
            ['SampleFieldState1'])
        self.assertEqual(comparison_groupings, {})

    def test_get_field_state_comparisons_no_comp_states(self):
        """get_field_state_comparisons() should raise a ValueError if no
        comparison field states are provided."""
        self.assertRaises(ValueError, get_field_state_comparisons,
                          self.dist_matrix_header, self.dist_matrix,
                          self.mapping_header, self.mapping, self.field,
                          [])

    def test_get_field_state_comparisons_bad_comp_state(self):
        """get_field_state_comparisons() should raise a ValueError if a
        non-existent comparison field state is provided."""
        self.assertRaises(ValueError, get_field_state_comparisons,
                          self.dist_matrix_header, self.dist_matrix,
                          self.mapping_header, self.mapping, self.field,
                          ['T0', 'Fast'])
        self.assertRaises(ValueError, get_field_state_comparisons,
                          self.dist_matrix_header, self.dist_matrix,
                          self.mapping_header, self.mapping, self.field,
                          ['Fast', 'T0'])

    def test_get_field_state_comparisons_invalid_distance_matrix(self):
        """Handles invalid distance matrix."""
        self.assertRaises(ValueError, get_field_state_comparisons,
                          ['Samp.1', 'Samp.2'],
                          array([[10.0, 0.0003], [0.0003, 0.0]]),
                          self.small_mapping_header, self.small_mapping,
                          self.small_field, ['SampleFieldState1'])

    def test_get_field_state_comparisons_no_shared_samples(self):
        """Handles invalid distance matrix."""
        self.assertRaises(ValueError, get_field_state_comparisons,
                          ['Samp.1', 'Samp.2'],
                          array([[10.0, 0.0003], [0.0003, 0.0]]),
                          self.small_mapping_header, [['foo', 'b', 'c'],
                          ['bar', 'bb', 'cc']],
                          self.small_field, ['SampleFieldState1'])

    def test_get_adjacent_distances(self):
        """ extracting adjacent distances works as expected
        """
        dm_str = ["\ts1\ts2\ts3", "s1\t0\t2\t4", "s2\t2\t0\t3.2",
                  "s3\t4\t3.2\t0"]
        dm_header, dm = parse_distmat(dm_str)
        # error cases: fewer than 2 valid sample ids
        self.assertRaises(ValueError,
                          get_adjacent_distances, dm_header, dm,
                          [])
        self.assertRaises(ValueError,
                          get_adjacent_distances, dm_header, dm,
                          ['s1'])
        self.assertRaises(ValueError,
                          get_adjacent_distances, dm_header, dm,
                          ['s0', 's1'])
        self.assertRaises(ValueError,
                          get_adjacent_distances, dm_header, dm,
                          ['s1', 's4'])

        # one pair of valid distances
        self.assertEqual(get_adjacent_distances(dm_header, dm, ['s1', 's2']),
                         ([2], [('s1', 's2')]))
        self.assertEqual(get_adjacent_distances(dm_header, dm, ['s1', 's1']),
                         ([0], [('s1', 's1')]))
        self.assertEqual(get_adjacent_distances(dm_header, dm, ['s1', 's3']),
                         ([4], [('s1', 's3')]))
        self.assertEqual(get_adjacent_distances(dm_header, dm, ['s2', 's3']),
                         ([3.2], [('s2', 's3')]))

        # multiple valid distances
        self.assertEqual(get_adjacent_distances(dm_header,
                                                dm,
                                                ['s1', 's2', 's3']),
                         ([2, 3.2], [('s1', 's2'), ('s2', 's3')]))
        self.assertEqual(get_adjacent_distances(dm_header,
                                                dm,
                                                ['s1', 's3', 's2', 's1']),
                         ([4, 3.2, 2], [('s1', 's3'), ('s3', 's2'), ('s2', 's1')]))

        # mixed valid and invalid distances ignores invalid distances
        self.assertEqual(get_adjacent_distances(dm_header,
                                                dm,
                                                ['s1', 's3', 's4', 's5', 's6', 's2', 's1']),
                         ([4, 3.2, 2], [('s1', 's3'), ('s3', 's2'), ('s2', 's1')]))
        # strict=True results in missing sample ids raising an error
        self.assertRaises(ValueError, get_adjacent_distances,
                          dm_header,
                          dm,
                          ['s1',
                           's3',
                           's4',
                           's5',
                           's6',
                           's2',
                           's1'],
                          strict=True)

    def test_get_ordered_coordinates(self):
        """get_ordered_coordinates functions as expected """
        pc_lines = ["Eigvals\t4",
                    "191.54\t169.99\t30.45\t19.19",
                    "",
                    "Proportion explained\t4",
                    "18.13\t16.09\t2.88\t1.66",
                    "",
                    "Species\t0\t0",
                    "",
                    "Site\t5\t4",
                    "s1\t-0.049\t0.245\t0.146\t-0.036",
                    "s5\t-0.267\t-0.228\t-0.024\t-0.095",
                    "s3\t-0.285\t-0.260\t-0.017\t-0.070",
                    "s2\t-0.002\t0.216\t-0.052\t-0.085",
                    "s4\t-0.328\t-0.299\t-0.025\t0.051",
                    "",
                    "Biplot\t0\t0",
                    "",
                    "Site constraints\t0\t0",
                    ""]

        pc = parse_coords(StringIO('\n'.join(pc_lines)))
        expected_coords = [[-0.049, 0.245, 0.146, -0.036],
                           [-0.002, 0.216, -0.052, -0.085],
                           [-0.285, -0.260, -0.017, -0.070],
                           [-0.328, -0.299, -0.025, 0.051],
                           [-0.267, -0.228, -0.024, -0.095]]
        expected_sids = ['s1', 's2', 's3', 's4', 's5']
        actual_coords, actual_sids = get_ordered_coordinates(
            pc[0], pc[1], ['s1', 's2', 's3', 's4', 's5'])
        assert_almost_equal(actual_coords, expected_coords)
        self.assertEqual(actual_sids, expected_sids)

        pc = parse_coords(StringIO('\n'.join(pc_lines)))
        expected_coords = [[-0.049, 0.245, 0.146, -0.036],
                           [-0.267, -0.228, -0.024, -0.095]]
        expected_sids = ['s1', 's5']
        actual_coords, actual_sids = get_ordered_coordinates(
            pc[0], pc[1], ['s1', 's5'])
        assert_almost_equal(actual_coords, expected_coords)
        self.assertEqual(actual_sids, expected_sids)

        pc = parse_coords(StringIO('\n'.join(pc_lines)))
        expected_coords = [[-0.049, 0.245, 0.146, -0.036],
                           [-0.267, -0.228, -0.024, -0.095]]
        expected_sids = ['s1', 's5']
        actual_coords, actual_sids = get_ordered_coordinates(
            pc[0], pc[1], ['s1', 's6', 's5'])
        assert_almost_equal(actual_coords, expected_coords)
        self.assertEqual(actual_sids, expected_sids)

        pc = parse_coords(StringIO('\n'.join(pc_lines)))
        expected_coords = [[-0.049, 0.245, 0.146, -0.036],
                           [-0.267, -0.228, -0.024, -0.095]]
        expected_sids = ['s1', 's5']
        self.assertRaises(ValueError, get_ordered_coordinates,
                          pc[0], pc[1], ['s1', 's6', 's5'], strict=True)

    def test_validate_input_bad_input(self):
        """_validate_input() should raise ValueErrors on bad input."""
        self.assertRaises(ValueError, _validate_input,
                          None, None, None, None, None)
        self.assertRaises(ValueError, _validate_input,
                          self.dist_matrix_header, self.dist_matrix,
                          self.mapping_header, self.mapping, None)
        self.assertRaises(ValueError, _validate_input,
                          self.dist_matrix_header, 12,
                          self.mapping_header, self.mapping, None)
        self.assertRaises(ValueError, _validate_input,
                          self.dist_matrix_header, self.dist_matrix,
                          self.mapping_header, self.mapping, 42)
        self.assertRaises(ValueError, _validate_input,
                          self.dist_matrix_header, self.dist_matrix,
                          self.mapping_header, self.mapping, "aeiou")

    def test_validate_input_good_input(self):
        """_validate_input() should not raise any errors on good input."""
        _validate_input(self.dist_matrix_header, self.dist_matrix,
                        self.mapping_header, self.mapping, "Treatment")

    def test_get_indices_several_existing_items(self):
        """_get_indices() should return a list of valid indices for several
        existing items."""
        control_ids = ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593']
        exp_control_indices = [0, 1, 2, 3, 4]

        fast_ids = ['PC.607', 'PC.634', 'PC.635', 'PC.636']
        exp_fast_indices = [5, 6, 7, 8]

        obs_control = _get_indices(self.dist_matrix_header, control_ids)
        self.assertEqual(obs_control, exp_control_indices)

        obs_fast = _get_indices(self.dist_matrix_header, fast_ids)
        self.assertEqual(obs_fast, exp_fast_indices)

    def test_get_indices_one_existing_item_list(self):
        """_get_indices() should return a list of size 1 for a single item in a
        list that exists in the search list."""
        item_to_find = ['PC.355']
        self.assertEqual(_get_indices(self.dist_matrix_header, item_to_find),
                         [1])

    def test_get_indices_one_existing_item_scalar(self):
        """_get_indices() should return a list of size 1 for a single item that
        exists in the search list."""
        item_to_find = 'PC.355'
        self.assertEqual(_get_indices(self.dist_matrix_header, item_to_find),
                         [1])

    def test_get_indices_no_existing_item(self):
        """_get_indices() should return an empty list if no items exist in the
        search list."""
        item_to_find = 'PC.4242'
        self.assertEqual(_get_indices(self.dist_matrix_header, item_to_find),
                         [])
        item_to_find = 42
        self.assertEqual(_get_indices(self.dist_matrix_header, item_to_find),
                         [])
        item_to_find = ['PC.4242', 'CP.2424']
        self.assertEqual(_get_indices(self.dist_matrix_header, item_to_find),
                         [])

        item_to_find = ['PC.4242', 'CP.2424', 56]
        self.assertEqual(_get_indices(self.dist_matrix_header, item_to_find),
                         [])

    def test_get_indices_no_items_to_search(self):
        """_get_indices() should return an empty list if no search items are
        given."""
        item_to_find = []
        self.assertEqual(_get_indices(self.dist_matrix_header, item_to_find),
                         [])
        item_to_find = ''
        self.assertEqual(_get_indices(self.dist_matrix_header, item_to_find),
                         [])
        item_to_find = None
        self.assertEqual(_get_indices(self.dist_matrix_header, item_to_find),
                         [])

    def test_get_indices_null_or_empty_search_list(self):
        """_get_indices() should throw an error if the search list is None, and
        return an empty list if the search list is empty."""
        search_list = None
        self.assertRaises(ValueError, _get_indices, search_list, 'item')

        search_list = []
        self.assertEqual(_get_indices(search_list, 'item'), [])

        search_list = ''
        self.assertEqual(_get_indices(search_list, 'item'), [])

    def test_get_groupings_no_field_states(self):
        """_get_groupings() should return an empty list if there are no field
        states in the groupings dictionary."""
        self.assertEqual(_get_groupings(self.dist_matrix_header,
                                        self.dist_matrix, {}, within=True), [])

        self.assertEqual(_get_groupings(self.dist_matrix_header,
                                        self.dist_matrix, {}, within=False), [])

    def test_get_groupings_within_tiny_dataset(self):
        """_get_groupings() should return an empty list for a single-sample
        dataset as the diagonal is omitted for within distances."""
        self.assertEqual(_get_groupings(self.tiny_dist_matrix_header,
                                        self.tiny_dist_matrix, self.tiny_groups, within=True), [])

    def test_get_groupings_between_tiny_dataset(self):
        """_get_groupings() should return an empty list for a single-sample
        dataset as there is only one field state, so no between distances can
        be computed."""
        self.assertEqual(_get_groupings(self.tiny_dist_matrix_header,
                                        self.tiny_dist_matrix, self.tiny_groups, within=False), [])

    def test_get_groupings_invalid_distance_matrix(self):
        """Handles asymmetric and/or hollow distance matrices correctly."""
        self.assertRaises(ValueError, _get_groupings, ['foo', 'bar'],
                          matrix([[0.0, 0.7], [0.7, 0.01]]), self.tiny_groups)

        # Should not raise error if we suppress the check.
        _get_groupings(['foo', 'bar'], matrix([[0.0, 0.7], [0.7, 0.01]]),
                       self.tiny_groups,
                       suppress_symmetry_and_hollowness_check=True)

    def test_extract_per_individual_states_from_sample_metadata(self):
        """extract_per_individual_states_from_sample_metadata functions as expected
        """
        expected = {'001': ['001A', '001B'],
                    '006': ['006A', '006B'],
                    '007': ['007A', '007B'],
                    '008': ['008A', '008B']}
        actual = extract_per_individual_states_from_sample_metadata(
            self.individual_states_and_responses_map_f1,
            state_category="TreatmentState",
            state_values=["Pre", "Post"],
            individual_identifier_category="PersonalID")
        self.assertEqual(actual, expected)

        # change to state_values order is reflected in order
        # of output sample ids
        expected = {'001': ['001B', '001A'],
                    '006': ['006B', '006A'],
                    '007': ['007B', '007A'],
                    '008': ['008B', '008A']}
        actual = extract_per_individual_states_from_sample_metadata(
            self.individual_states_and_responses_map_f1,
            state_category="TreatmentState",
            state_values=["Post", "Pre"],
            individual_identifier_category="PersonalID")
        self.assertEqual(actual, expected)

        # don't filter missing data
        expected = {'001': ['001A', '001B'],
                    '006': ['006A', '006B'],
                    '007': ['007A', '007B'],
                    '008': ['008A', '008B'],
                    '009': [None, 'post.only'],
                    '010': ['pre.only', None]}
        actual = extract_per_individual_states_from_sample_metadata(
            self.individual_states_and_responses_map_f1,
            state_category="TreatmentState",
            state_values=["Pre", "Post"],
            individual_identifier_category="PersonalID",
            filter_missing_data=False)
        self.assertEqual(actual, expected)

        # alt input file with more states
        expected = {'001': ['001A', '001B', '001C']}
        actual = extract_per_individual_states_from_sample_metadata(
            self.individual_states_and_responses_map_f2,
            state_category="TreatmentState",
            state_values=["Pre", "Post", "PostPost"],
            individual_identifier_category="PersonalID")
        self.assertEqual(actual, expected)

        # unlisted states are ignored
        expected = {'001': ['001A', '001B']}
        actual = extract_per_individual_states_from_sample_metadata(
            self.individual_states_and_responses_map_f2,
            state_category="TreatmentState",
            state_values=["Pre", "Post"],
            individual_identifier_category="PersonalID")
        self.assertEqual(actual, expected)

    def test_extract_per_individual_states_from_sample_metadata_invalid(self):
        """extract_per_individual_states_from_sample_metadata handles invalid input
        """
        self.assertRaises(KeyError,
                          extract_per_individual_states_from_sample_metadata,
                          self.individual_states_and_responses_map_f1,
                          state_category="some-invalid-category",
                          state_values=["Pre", "Post"],
                          individual_identifier_category="PersonalID")
        self.assertRaises(KeyError,
                          extract_per_individual_states_from_sample_metadata,
                          self.individual_states_and_responses_map_f1,
                          state_category="TreatmentState",
                          state_values=["Pre", "Post"],
                          individual_identifier_category="some-other-invalid-category")

    def test_extract_per_individual_state_metadatum_from_sample_metadata(self):
        """extract_per_individual_state_metadatum_from_sample_metadata functions as expected
        """
        veil_expected = {'001': [6.9, 9.3],
                         '006': [4.2, 5.1],
                         '007': [12.0, 1.8],
                         '008': [10.0, None]}
        actual = extract_per_individual_state_metadatum_from_sample_metadata(
            self.individual_states_and_responses_map_f1,
            state_category="TreatmentState",
            state_values=["Pre", "Post"],
            individual_identifier_category="PersonalID",
            metadata_category="VeillonellaAbundance",
            process_f=float)
        self.assertEqual(actual, veil_expected)

        # different metadata_category yields different result
        strep_expected = {'001': [57.4, 26],
                          '006': [19, 15.2],
                          '007': [33.2, 50],
                          '008': [3.2, 20]}
        actual = extract_per_individual_state_metadatum_from_sample_metadata(
            self.individual_states_and_responses_map_f1,
            state_category="TreatmentState",
            state_values=["Pre", "Post"],
            individual_identifier_category="PersonalID",
            metadata_category="StreptococcusAbundance",
            process_f=float)
        self.assertEqual(actual, strep_expected)

        # different metadata_category yields different result
        response_expected = {'001': ["Improved", "Improved"],
                             '006': ["Improved", "Improved"],
                             '007': ["Worsened", "Worsened"],
                             '008': ["Worsened", "Worsened"]}
        actual = extract_per_individual_state_metadatum_from_sample_metadata(
            self.individual_states_and_responses_map_f1,
            state_category="TreatmentState",
            state_values=["Pre", "Post"],
            individual_identifier_category="PersonalID",
            metadata_category="Response",
            process_f=str)
        self.assertEqual(actual, response_expected)

        # alt input file with more states
        expected = {'001': [6.9, 9.3, 10.1]}
        actual = extract_per_individual_state_metadatum_from_sample_metadata(
            self.individual_states_and_responses_map_f2,
            state_category="TreatmentState",
            state_values=["Pre", "Post", "PostPost"],
            individual_identifier_category="PersonalID",
            metadata_category="VeillonellaAbundance",
            process_f=float)
        self.assertEqual(actual, expected)

        # unlisted states are ignored
        expected = {'001': [6.9, 9.3]}
        actual = extract_per_individual_state_metadatum_from_sample_metadata(
            self.individual_states_and_responses_map_f2,
            state_category="TreatmentState",
            state_values=["Pre", "Post"],
            individual_identifier_category="PersonalID",
            metadata_category="VeillonellaAbundance",
            process_f=float)
        self.assertEqual(actual, expected)

    def test_extract_per_individual_state_metadatum_from_sample_metadata_invalid(
            self):
        """extract_per_individual_state_metadatum_from_sample_metadata handles invalid column header
        """
        self.assertRaises(KeyError,
                          extract_per_individual_state_metadatum_from_sample_metadata,
                          self.individual_states_and_responses_map_f1,
                          state_category="TreatmentState",
                          state_values=["Pre", "Post"],
                          individual_identifier_category="PersonalID",
                          metadata_category="some-non-existant-category",
                          process_f=float)

    def test_extract_per_individual_state_metadata_from_sample_metadata_and_biom(
            self):
        """extract_per_individual_state_metadata_from_sample_metadata_and_biom functions as expected
        """
        # single observations
        o1_expected = {'o1': {'001': [22, 10],
                              '006': [25, 4],
                              '007': [33, 26],
                              '008': [99, 75]}}
        actual = extract_per_individual_state_metadata_from_sample_metadata_and_biom(
            self.individual_states_and_responses_map_f1,
            self.paired_difference_biom1,
            state_category="TreatmentState",
            state_values=["Pre", "Post"],
            individual_identifier_category="PersonalID",
            observation_ids=['o1'])
        self.assertEqual(actual, o1_expected)

        # all observations
        all_expected = {'o1': {'001': [22, 10],
                               '006': [25, 4],
                               '007': [33, 26],
                               '008': [99, 75]},
                        'o2': {'001': [10, 44],
                               '006': [3, 99],
                               '007': [8, 18],
                               '008': [64, 164]},
                        'o3': {'001': [10, 50],
                               '006': [50, 10],
                               '007': [10, 50],
                               '008': [50, 10]}}
        actual = extract_per_individual_state_metadata_from_sample_metadata_and_biom(
            self.individual_states_and_responses_map_f1,
            self.paired_difference_biom1,
            state_category="TreatmentState",
            state_values=["Pre", "Post"],
            individual_identifier_category="PersonalID",
            observation_ids=None)
        self.assertEqual(actual, all_expected)

        # invalid observation id
        self.assertRaises(
            UnknownIDError,
            extract_per_individual_state_metadata_from_sample_metadata_and_biom,
            self.individual_states_and_responses_map_f1,
            self.paired_difference_biom1,
            state_category="TreatmentState",
            state_values=["Pre", "Post"],
            individual_identifier_category="PersonalID",
            observation_ids=['o1', 'bad.obs.id'])

    def test_group_by_sample_metadata(self):
        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        collapsed_md = _collapse_metadata(in_f, ['replicate-group', 'subject'])
        actual = _group_by_sample_metadata(collapsed_md)
        expected1 = {('1', '1'): set(('f1', 'f2')),
                     ('2', '1'): set(('f5', 'f6', 'p1')),
                     ('3', '1'): set(('not16S.1', )),
                     ('1', '2'): set(('f3', 'f4')),
                     ('2', '2'): set(('p2', 't1', 't2'))}
        expected2 = {'f1': ('1', '1'),
                     'f2': ('1', '1'),
                     'f5': ('2', '1'),
                     'f6': ('2', '1'),
                     'p1': ('2', '1'),
                     'not16S.1': ('3', '1'),
                     'f3': ('1', '2'),
                     'f4': ('1', '2'),
                     'p2': ('2', '2'),
                     't1': ('2', '2'),
                     't2': ('2', '2')}
        self.assertEqual(actual, (expected1, expected2))

        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        collapsed_md = _collapse_metadata(in_f, ['replicate-group'])
        actual = _group_by_sample_metadata(collapsed_md)
        expected1 = {('1', ): set(('f1', 'f2', 'f3', 'f4')),
                     ('2', ): set(('f5', 'f6', 'p1', 'p2', 't1', 't2')),
                     ('3', ): set(('not16S.1', ))}
        expected2 = {'f1': ('1', ), 'f2': ('1', ), 'f5': ('2', ), 'f6': ('2', ),
                     'p1': ('2', ), 'not16S.1': ('3', ), 'f3': ('1', ),
                     'f4': ('1', ), 'p2': ('2', ), 't1': ('2', ), 't2': ('2', )}
        self.assertEqual(actual, (expected1, expected2))

        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        collapsed_md = _collapse_metadata(in_f, ['subject'])
        actual = _group_by_sample_metadata(collapsed_md)
        expected1 = {('1', ): set(('f1', 'f2', 'f5', 'f6', 'p1', 'not16S.1')),
                     ('2', ): set(('f3', 'f4', 'p2', 't1', 't2'))}
        expected2 = {'f1': ('1', ), 'f2': ('1', ), 'f5': ('1', ), 'f6': ('1', ),
                     'p1': ('1', ), 'not16S.1': ('1', ), 'f3': ('2', ),
                     'f4': ('2', ), 'p2': ('2', ), 't1': ('2', ), 't2': ('2', )}
        self.assertEqual(actual, (expected1, expected2))

    def test_sample_id_from_group_id(self):
        sid_to_group_id1 = {'f1': (1, ), 'f2': (2, ), 'f5': (4, )}
        md = {}
        self.assertEqual(_sample_id_from_group_id('f1', md, sid_to_group_id1),
                         '1')
        self.assertEqual(_sample_id_from_group_id('f2', md, sid_to_group_id1),
                         '2')
        self.assertEqual(_sample_id_from_group_id('f5', md, sid_to_group_id1),
                         '4')

        sid_to_group_id2 = {'f1': (1, 1), 'f2': (1, 1), 'f5': (2, 1)}
        self.assertEqual(_sample_id_from_group_id('f1', md, sid_to_group_id2),
                         '1.1')
        self.assertEqual(_sample_id_from_group_id('f2', md, sid_to_group_id2),
                         '1.1')
        self.assertEqual(_sample_id_from_group_id('f5', md, sid_to_group_id2),
                         '2.1')

        sid_to_group_id3 = {'f1': (1, 1, 2), 'f5': (2, 1, 0)}
        self.assertEqual(_sample_id_from_group_id('f1', md, sid_to_group_id3),
                         '1.1.2')
        self.assertEqual(_sample_id_from_group_id('f5', md, sid_to_group_id3),
                         '2.1.0')

        self.assertRaises(KeyError, _sample_id_from_group_id, 'f2', md,
                          sid_to_group_id3)

    def test_collapse_samples(self):
        """Collapsing samples functions as expected
        """
        # #OTU ID	f1	f2	f3
        # o1	0.0	1.0	2.0
        # o2	3.0	4.0	5.0
        t1 = Table(np.array([[0, 1, 2], [3, 4, 5]]),
                   ['o1', 'o2'], ['f1', 'f2', 'f3'])
        collapse_fields = ['replicate-group', 'subject']
        for e in get_collapse_fns().keys():
            # all collapse functions work without failure
            in_f = StringIO(self._group_by_sample_metadata_map_f1)
            collapse_samples(t1, in_f, collapse_fields, e)
        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        self.assertRaises(KeyError, collapse_samples, t1, in_f,
                          collapse_fields, "not-a-valid-mode")

        # test with a few collapse functions (the collapse functions
        # are tested individually, so don't need to test them all here)
        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        md, t = collapse_samples(t1, in_f, collapse_fields, 'sum')
        self.assertEqual(t.get_value_by_ids('o1', '1.1'), 1.0)
        self.assertEqual(t.get_value_by_ids('o1', '1.2'), 2.0)
        self.assertEqual(t.get_value_by_ids('o2', '1.1'), 7.0)
        self.assertEqual(t.get_value_by_ids('o2', '1.2'), 5.0)

        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        md, t = collapse_samples(t1, in_f, collapse_fields, 'mean')
        self.assertEqual(t.get_value_by_ids('o1', '1.1'), 0.5)
        self.assertEqual(t.get_value_by_ids('o1', '1.2'), 2.0)
        self.assertEqual(t.get_value_by_ids('o2', '1.1'), 3.5)
        self.assertEqual(t.get_value_by_ids('o2', '1.2'), 5.0)

        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        md, t = collapse_samples(t1, in_f, collapse_fields, 'first')
        self.assertEqual(t.get_value_by_ids('o1', '1.1'), 0.0)
        self.assertEqual(t.get_value_by_ids('o1', '1.2'), 2.0)
        self.assertEqual(t.get_value_by_ids('o2', '1.1'), 3.0)
        self.assertEqual(t.get_value_by_ids('o2', '1.2'), 5.0)

    def test_collapse_to_first(self):
        """ Table collapse function _collapse_to_first functions as expected
        """
        # #OTU ID	s1	s2	s3
        # o1	0.0	1.0	2.0
        # o2	3.0	4.0	5.0
        t1 = Table(np.array([[0, 1, 2], [3, 4, 5]]),
                   ['o1', 'o2'], ['s1', 's2', 's3'])
        self.assertEqual(list(_collapse_to_first(t1, "observation")),
                         [0.0, 3.0])
        self.assertEqual(list(_collapse_to_first(t1, "sample")),
                         [0.0, 1.0, 2.0])

    def test_collapse_to_median(self):
        """ Table collapse function _collapse_to_median functions as expected
        """
        # #OTU ID	s1	s2	s3
        # o1	0.0	1.0	2.0
        # o2	3.0	4.0	5.0
        t1 = Table(np.array([[0, 1, 2], [3, 4, 5]]),
                   ['o1', 'o2'], ['s1', 's2', 's3'])
        self.assertEqual(list(_collapse_to_median(t1, "observation")),
                         [1.0, 4.0])
        self.assertEqual(list(_collapse_to_median(t1, "sample")),
                         [1.5, 2.5, 3.5])

    def test_collapse_to_random(self):
        """ Table collapse function _collapse_to_random functions as expected
        """
        # #OTU ID	s1	s2	s3
        # o1	0.0	1.0	2.0
        # o2	3.0	4.0	5.0
        t1 = Table(np.array([[0, 1, 2], [3, 4, 5]]),
                   ['o1', 'o2'], ['s1', 's2', 's3'])
        e = list(_collapse_to_random(t1, "observation"))
        self.assertTrue(e[0] in [0, 1, 2])
        self.assertTrue(e[1] in [3, 4, 5])

        e = list(_collapse_to_random(t1, "sample"))
        self.assertTrue(e[0] in [0, 3])
        self.assertTrue(e[1] in [1, 4])
        self.assertTrue(e[2] in [2, 5])

    def test_collapse_to_sum(self):
        """ Table collapse function _collapse_to_sum functions as expected
        """
        # #OTU ID	s1	s2	s3
        # o1	0.0	1.0	2.0
        # o2	3.0	4.0	5.0
        t1 = Table(np.array([[0, 1, 2], [3, 4, 5]]),
                   ['o1', 'o2'], ['s1', 's2', 's3'])
        self.assertEqual(list(_collapse_to_sum(t1, "observation")),
                         [3.0, 12.0])
        self.assertEqual(list(_collapse_to_sum(t1, "sample")),
                         [3.0, 5.0, 7.0])

    def test_collapse_to_mean(self):
        """ Table collapse function _collapse_to_mean functions as expected
        """
        # #OTU ID	s1	s2	s3
        # o1	0.0	1.0	2.0
        # o2	3.0	4.0	5.0
        t1 = Table(np.array([[0, 1, 2], [3, 4, 5]]),
                   ['o1', 'o2'], ['s1', 's2', 's3'])
        self.assertEqual(list(_collapse_to_mean(t1, "observation")),
                         [1.0, 4.0])
        self.assertEqual(list(_collapse_to_mean(t1, "sample")),
                         [1.5, 2.5, 3.5])

    def test_collapse_metadata(self):
        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        actual = _collapse_metadata(
            in_f, ['replicate-group', 'subject'])
        # correct collapsing
        self.assertEqual(actual['SampleID'][('1', '1')], ('f1', 'f2'))
        self.assertEqual(actual['SampleID'][('1', '2')], ('f3', 'f4'))
        self.assertEqual(actual['SampleID'][('2', '1')], ('f5', 'f6', 'p1'))
        self.assertEqual(actual['SampleID'][('2', '2')], ('p2', 't1', 't2'))
        self.assertEqual(actual['SampleID'][('3', '1')], ('not16S.1', ))
        # original values tuple-ized
        self.assertEqual(actual['BarcodeSequence'][('1', '1')],
                         ('ACACTGTTCATG', 'ACCAGACGATGC'))
        self.assertEqual(actual['BarcodeSequence'][('1', '2')],
                         ('ACCAGACGATGC', 'ACCAGACGATGC'))

        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        self.assertRaises(KeyError, _collapse_metadata, in_f,
                          ['not-a-header'])

    def test_mapping_lines_from_collapsed_df(self):
        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        collapsed_df = _collapse_metadata(
        in_f, ['subject'])
        expected = mapping_lines_from_collapsed_df(collapsed_df)
        self.assertTrue(expected[0].startswith("#SampleID	original-sample-ids	BarcodeSequence	LinkerPrimerSequence"))
        self.assertTrue(expected[1].startswith('1	(f1, f2, f5, f6, p1, not16S.1)	'))
        self.assertTrue(expected[2].startswith('2	(f3, f4, p2, t1, t2)	'))

        in_f = StringIO(self._group_by_sample_metadata_map_f1)
        collapsed_df = _collapse_metadata(
            in_f, ['replicate-group', 'subject'])
        expected = mapping_lines_from_collapsed_df(collapsed_df)
        self.assertTrue(expected[0].startswith("#SampleID	original-sample-ids	BarcodeSequence	LinkerPrimerSequence"))
        self.assertTrue(expected[1].startswith('1.1	(f1, f2)	(ACACTGTTCATG, ACCAGACGATGC)	GTGCCAGCMGCCGCGGTAA	feces'))
        self.assertTrue(expected[2].startswith('1.2	(f3, f4)	ACCAGACGATGC	GTGCCAGCMGCCGCGGTAA	feces'))


individual_states_and_responses_map_f1 = """#SampleID	PersonalID	Response	TreatmentState	StreptococcusAbundance	VeillonellaAbundance
001A	001	Improved	Pre	57.4	6.9
001B	001	Improved	Post	26	9.3
006A	006	Improved	Pre	19	4.2
006B	006	Improved	Post	15.2	5.1
007A	007	Worsened	Pre	33.2	12
007B	007	Worsened	Post	50	1.8
008A	008	Worsened	Pre	3.2	10
008B	008	Worsened	Post	20	n/a
post.only	009	Worsened	Post	22	42.0
pre.only	010	Worsened	Pre	21	41.0
"""

paired_difference_biom_f1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-07-10T09:22:39.602392","matrix_type": "sparse","matrix_element_type": "float","shape": [3, 10],"data": [[0,0,22.0],[0,1,10.0],[0,2,25.0],[0,3,4.0],[0,4,33.0],[0,5,26.0],[0,6,99.0],[0,7,75.0],[0,8,66.0],[0,9,67.0],[1,0,10.0],[1,1,44.0],[1,2,3.0],[1,3,99.0],[1,4,8.0],[1,5,18.0],[1,6,64.0],[1,7,164.0],[1,8,22.0],[1,9,22.0],[2,0,10.0],[2,1,50.0],[2,2,50.0],[2,3,10.0],[2,4,10.0],[2,5,50.0],[2,6,50.0],[2,7,10.0],[2,8,10.0],[2,9,50.0]],"rows": [{"id": "o1", "metadata": null},{"id": "o2", "metadata": null},{"id": "o3", "metadata": null}],"columns": [{"id": "001A", "metadata": null},{"id": "001B", "metadata": null},{"id": "006A", "metadata": null},{"id": "006B", "metadata": null},{"id": "007A", "metadata": null},{"id": "007B", "metadata": null},{"id": "008A", "metadata": null},{"id": "008B", "metadata": null},{"id": "post.only", "metadata": null},{"id": "pre.only", "metadata": null}]}"""

individual_states_and_responses_map_f2 = """#SampleID	PersonalID	Response	TreatmentState	StreptococcusAbundance	VeillonellaAbundance
001A	001	Improved	Pre	57.4	6.9
001B	001	Improved	Post	26	9.3
001C	001	Improved	PostPost	22	10.1
"""

_group_by_sample_metadata_map_f1 = """#SampleID	BarcodeSequence	LinkerPrimerSequence	SampleType	year	month	day	subject	replicate-group	days_since_epoch	Description
f1	ACACTGTTCATG	GTGCCAGCMGCCGCGGTAA	feces	2008	10	22	1	1	14174	fecal1
f2	ACCAGACGATGC	GTGCCAGCMGCCGCGGTAA	feces	2008	10	23	1	1	14175	fecal2
f3	ACCAGACGATGC	GTGCCAGCMGCCGCGGTAA	feces	2008	10	23	2	1	14175	identical sequences to fecal2
f4	ACCAGACGATGC	GTGCCAGCMGCCGCGGTAA	feces	2008	10	23	2	1	14175	all sequences identical, map to GG 295053 at 97 percent id
f5	ACCAGACGATGC	GTGCCAGCMGCCGCGGTAA	feces	2008	10	23	1	2	14175	derived from f3 with some changes to sequences to add one new otu
f6	ACCAGACGATGC	GTGCCAGCMGCCGCGGTAA	feces	2008	10	23	1	2	14175	derived from f4 with some changes to sequences to add one new otu
p1	AACGCACGCTAG	GTGCCAGCMGCCGCGGTAA	L_palm	2008	10	21	1	2	14173	palm1, contains one randomly generated sequence
p2	ACACTGTTCATG	GTGCCAGCMGCCGCGGTAA	L_palm	2008	10	22	2	2	14174	palm2
t1	AGTGAGAGAAGC	GTGCCAGCMGCCGCGGTAA	Tongue	2008	10	21	2	2	14173	tongue1, contains one randomly generated sequence
t2	ATACTATTGCGC	GTGCCAGCMGCCGCGGTAA	Tongue	2008	10	22	2	2	14174	tongue2
not16S.1	ATACTATTGCGC	GTGCCAGCMGCCGCGGTAA	Other	2008	10	22	1	3	14174	randomly generated sequence plus some variants, these should not map to 16S
"""

if __name__ == '__main__':
    main()
