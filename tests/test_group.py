#!/usr/bin/env python

"""Tests public and private functions in the group module."""

__author__ = "Jai Rideout"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jai Rideout", "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Jai Rideout"
__email__ = "jr378@nau.edu"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_mapping_file, parse_distmat, group_by_field
from qiime.group import get_grouped_distances, get_all_grouped_distances,\
    _get_indices, _get_groupings, _validate_input

class GroupTests(TestCase):
    """Tests of the group module."""
    def setUp(self):
        """Create some data to be used in the tests."""
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

        # Parse mapping "file" (faked here).
        self.mapping, self.mapping_header, self.comments = parse_mapping_file(
                self.mapping_string)
        mapping_data = [self.mapping_header]
        mapping_data.extend(self.mapping)
        self.groups = group_by_field(mapping_data, self.field)

        # Parse distance matrix "file" (faked here).
        self.dist_matrix_header, self.dist_matrix = parse_distmat(
                self.dist_matrix_string)

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

if __name__ == '__main__':
    main()
