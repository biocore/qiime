#!/usr/bin/env python

"""Tests public and private functions in the group module."""

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jai Ram Rideout",
               "Greg Caporaso",
               "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Release"

from numpy import array, matrix
from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_mapping_file, parse_distmat, group_by_field, parse_coords
from qiime.group import (get_grouped_distances, get_all_grouped_distances,
    get_field_state_comparisons, _get_indices, _get_groupings, _validate_input,
    get_adjacent_distances, get_ordered_coordinates)

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

    def test_get_grouped_distances_within(self):
        """get_grouped_distances() should return a list of within distance
        groupings."""
        groupings = get_grouped_distances(self.dist_matrix_header,
            self.dist_matrix, self.mapping_header, self.mapping,
            self.field, within=True)
        expected = [
            ('Control', 'Control', [0.625, 0.623, 0.60999999999999999, \
                                    0.57699999999999996, 0.61499999999999999, \
                                    0.64200000000000002, 0.67300000000000004, \
                                    0.68200000000000005, 0.73699999999999999, \
                                    0.70399999999999996]),
            ('Fast', 'Fast', [0.71799999999999997, 0.66600000000000004, \
                              0.72699999999999998, 0.59999999999999998, \
                              0.57799999999999996, 0.623])]
        self.assertEqual(groupings, expected)

    def test_get_grouped_distances_between(self):
        """get_grouped_distances() should return a list of between distance
        groupings."""
        groupings = get_grouped_distances(self.dist_matrix_header,
            self.dist_matrix, self.mapping_header, self.mapping,
            self.field, within=False)
        expected = [
            ('Control', 'Fast', [0.72899999999999998, 0.80000000000000004, \
                                 0.72099999999999997, 0.76500000000000001, \
                                 0.77600000000000002, 0.74399999999999999, \
                                 0.749, 0.67700000000000005, \
                                 0.73399999999999999, 0.77700000000000002, \
                                 0.73299999999999998, 0.72399999999999998, \
                                 0.69599999999999995, 0.67500000000000004, \
                                 0.65400000000000003, 0.69599999999999995, \
                                 0.73099999999999998, 0.75800000000000001, \
                                 0.73799999999999999, 0.73699999999999999])]
        self.assertEqual(groupings, expected)

    def test_get_all_grouped_distances_within(self):
        """get_all_grouped_distances() should return a list of distances for
        all samples with the same field value."""
        groupings = get_all_grouped_distances(self.dist_matrix_header,
            self.dist_matrix, self.mapping_header, self.mapping,
            self.field, within=True)
        expected =  [0.625, 0.623, 0.60999999999999999, 0.57699999999999996,
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
        self.assertFloatEqual(comparison_groupings, expected)

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
        self.assertFloatEqual(comparison_groupings, expected)

    def test_get_field_state_comparisons_small(self):
        """get_field_state_comparisons() should return a 2D dictionary of
        distances between a field state and its comparison field states."""
        comparison_groupings = get_field_state_comparisons(
                self.small_dist_matrix_header, self.small_dist_matrix,
                self.small_mapping_header, self.small_mapping,
                self.small_field, ['SampleFieldState1'])
        expected = {'SampleFieldState2': {'SampleFieldState1': [0.5]}}
        self.assertFloatEqual(comparison_groupings, expected)

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

    def test_get_adjacent_distances(self):
        """ extracting adjacent distances works as expected
        """
        dm_str = ["\ts1\ts2\ts3", "s1\t0\t2\t4", "s2\t2\t0\t3.2",
                        "s3\t4\t3.2\t0"]
        dm_header, dm = parse_distmat(dm_str)
        # error cases: fewer than 2 valid sample ids
        self.assertRaises(ValueError,
                          get_adjacent_distances,dm_header, dm,
                          [])
        self.assertRaises(ValueError,
                          get_adjacent_distances,dm_header, dm,
                          ['s1'])
        self.assertRaises(ValueError,
                          get_adjacent_distances,dm_header, dm,
                          ['s0','s1'])
        self.assertRaises(ValueError,
                          get_adjacent_distances,dm_header, dm,
                          ['s1','s4'])
        
        # one pair of valid distances
        self.assertEqual(get_adjacent_distances(dm_header, dm, ['s1','s2']),
                         ([2],[('s1','s2')]))
        self.assertEqual(get_adjacent_distances(dm_header, dm, ['s1','s1']),
                         ([0],[('s1','s1')]))
        self.assertEqual(get_adjacent_distances(dm_header, dm, ['s1','s3']),
                         ([4],[('s1','s3')]))
        self.assertEqual(get_adjacent_distances(dm_header, dm, ['s2','s3']),
                         ([3.2],[('s2','s3')]))
        
        # multiple valid distances
        self.assertEqual(get_adjacent_distances(dm_header, 
                                                dm, 
                                                ['s1','s2','s3']),
                         ([2,3.2],[('s1','s2'),('s2','s3')]))
        self.assertEqual(get_adjacent_distances(dm_header, 
                                                dm, 
                                                ['s1','s3','s2','s1']),
                         ([4,3.2,2],[('s1','s3'),('s3','s2'),('s2','s1')]))
        
        # mixed valid and invalid distances ignores invalid distances
        self.assertEqual(get_adjacent_distances(dm_header, 
                                                dm, 
                                                ['s1','s3','s4','s5','s6','s2','s1']),
                         ([4,3.2,2],[('s1','s3'),('s3','s2'),('s2','s1')]))
        # strict=True results in missing sample ids raising an error
        self.assertRaises(ValueError,get_adjacent_distances,
                                     dm_header, 
                                     dm,
                                     ['s1','s3','s4','s5','s6','s2','s1'],
                                     strict=True)
                                     
    def test_get_ordered_coordinates(self):
        """get_ordered_coordinates functions as expected """
        pc_lines = ["pc vector number\t1\t2\t3\t4",
                    "s1\t-0.049\t0.245\t0.146\t-0.036",
                    "s5\t-0.267\t-0.228\t-0.024\t-0.095",
                    "s3\t-0.285\t-0.260\t-0.017\t-0.070",
                    "s2\t-0.002\t0.216\t-0.052\t-0.085",
                    "s4\t-0.328\t-0.299\t-0.025\t0.051",
                    "",
                    "",
                    "eigvals\t191.54\t169.99\t30.45\t19.19",
                    "%% variation explained\t18.13\t16.09\t2.88\t1.66"]

        pc = parse_coords(pc_lines)
        expected_coords = [[-0.049, 0.245, 0.146, -0.036],
                           [-0.002, 0.216, -0.052, -0.085],
                           [-0.285, -0.260, -0.017, -0.070],
                           [-0.328, -0.299, -0.025, 0.051],
                           [-0.267, -0.228, -0.024, -0.095]]
        expected_sids = ['s1','s2','s3','s4','s5']
        actual_coords, actual_sids = get_ordered_coordinates(
         pc[0],pc[1],['s1','s2','s3','s4','s5'])
        self.assertEqual(actual_coords,expected_coords)
        self.assertEqual(actual_sids,expected_sids)

        pc = parse_coords(pc_lines)
        expected_coords = [[-0.049, 0.245, 0.146, -0.036],
                           [-0.267, -0.228, -0.024, -0.095]]
        expected_sids = ['s1','s5']
        actual_coords, actual_sids = get_ordered_coordinates(
         pc[0],pc[1],['s1','s5'])
        self.assertEqual(actual_coords,expected_coords)
        self.assertEqual(actual_sids,expected_sids)

        pc = parse_coords(pc_lines)
        expected_coords = [[-0.049, 0.245, 0.146, -0.036],
                           [-0.267, -0.228, -0.024, -0.095]]
        expected_sids = ['s1','s5']
        actual_coords, actual_sids = get_ordered_coordinates(
         pc[0],pc[1],['s1','s6','s5'])
        self.assertEqual(actual_coords,expected_coords)
        self.assertEqual(actual_sids,expected_sids)

        pc = parse_coords(pc_lines)
        expected_coords = [[-0.049, 0.245, 0.146, -0.036],
                           [-0.267, -0.228, -0.024, -0.095]]
        expected_sids = ['s1','s5']
        self.assertRaises(ValueError,get_ordered_coordinates,
                          pc[0],pc[1],['s1','s6','s5'],strict=True)
        
        

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
        exp_control_indices = [0,1,2,3,4]
        
        fast_ids = ['PC.607', 'PC.634', 'PC.635', 'PC.636']
        exp_fast_indices = [5,6,7,8]
        
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


if __name__ == '__main__':
    main()
