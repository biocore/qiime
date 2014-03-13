#!/usr/bin/env python
# File created on 21 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from numpy import array
from unittest import TestCase, main
from qiime.parse import parse_coords
from qiime.transform_coordinate_matrices import map_sample_ids, reorder_coords,\
    filter_coords_matrix, pad_coords_matrix, pad_coords_matrices,\
    get_procrustes_results, procrustes_monte_carlo


class ProcrustesTests(TestCase):

    """ Tests of the Procrustes wrapper code """

    def setUp(self):
        """ """
        self.pcoa1_f = pcoa1_f.split('\n')
        self.sample_ids1, self.coords1, self.eigvals1, self.pct_var1 =\
            parse_coords(self.pcoa1_f)
        self.pcoa2_f = pcoa2_f.split('\n')
        self.sample_ids2, self.coords2, self.eigvals2, self.pct_var2 =\
            parse_coords(self.pcoa2_f)
        self.pcoa3_f = pcoa3_f.split('\n')
        self.sample_ids3, self.coords3, self.eigvals3, self.pct_var3 =\
            parse_coords(self.pcoa3_f)
        self.pcoa4_f = pcoa4_f.split('\n')
        self.sample_ids4, self.coords4, self.eigvals4, self.pct_var4 =\
            parse_coords(self.pcoa3_f)

        self.sample_id_map1 = sample_id_map1

    def test_map_sample_ids(self):
        """Mapping and reordering of sample IDs functions as expected
        """
        expected = ['s4', 's1', 's2', 's3']
        actual = map_sample_ids(self.sample_ids1, self.sample_id_map1)
        self.assertEqual(actual, expected)
        expected = ['s1', 's2', 's3', 's4']
        actual = map_sample_ids(self.sample_ids2, self.sample_id_map1)
        self.assertEqual(actual, expected)
        expected = ['s1', 's2', 's3', 's4']
        actual = map_sample_ids(self.sample_ids3, self.sample_id_map1)
        self.assertEqual(actual, expected)
        # bad sample id raises KeyError
        self.assertRaises(KeyError, map_sample_ids, ['abcd', 'aaa', 'ccc'],
                          self.sample_id_map1)

    def test_reorder_coords(self):
        """Reordering of columns functions as expected """
        m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        in_sids = ['A', 'B', 'C']
        order = ['A', 'B', 'C']
        expected = m
        self.assertEqual(reorder_coords(m, in_sids, order), expected)

        in_sids = ['C', 'B', 'A']
        expected = [[7, 8, 9], [4, 5, 6], [1, 2, 3]]
        self.assertEqual(reorder_coords(m, in_sids, order), expected)

        in_sids = ['C', 'B', 'A']
        expected = [[7, 8, 9], [4, 5, 6], [1, 2, 3]]
        self.assertEqual(reorder_coords(m, in_sids, order), expected)

        # order leaves some samples out
        m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        in_sids = ['A', 'B', 'C']
        order = ['A', 'B']
        expected = [[1, 2, 3], [4, 5, 6]]
        self.assertEqual(reorder_coords(m, in_sids, order), expected)

    def test_reorder_coords_errors(self):
        """Reordering of columns handles errors """

        m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        in_sids = ['A', 'B', 'C']
        order = ['A', 'B', 'C', 'D']
        self.assertRaises(ValueError, reorder_coords, m, in_sids, order)

    def test_filter_coords_matrix(self):
        """ Removing coordinate from each sample functions as expected
        """
        m = array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        expected = array([[1], [4], [7], [10]])
        self.assertEqual(filter_coords_matrix(m, 1), expected)
        expected = array([[1, 2], [4, 5], [7, 8], [10, 11]])
        self.assertEqual(filter_coords_matrix(m, 2), expected)

    def test_pad_coords_matrix(self):
        """ padding a coordinates matrix with zeros functions as expected
        """
        m = array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        actual = pad_coords_matrix(m, 0)
        self.assertEqual(actual, m)

        m = array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        expected = array([[1, 2, 3, 0., 0.], [4, 5, 6, 0., 0.],
                          [7, 8, 9, 0., 0.], [10, 11, 12, 0., 0.]])
        actual = pad_coords_matrix(m, 2)
        self.assertEqual(actual, expected)

        m = array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        self.assertRaises(ValueError, pad_coords_matrix, m, -3)

    def test_pad_coords_matrices(self):
        """ padding the correct choice of a coords matrix with zeros functions
        """
        # first is longer
        m1 = array(
            [[1, 0, 3, 9, 4], [4, 5, 6, 94, 5], [7, 8, 9, 6, 6], [10, 11, 12, 4, 6]])
        m2 = array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        m2_expected = array([[1, 2, 3, 0., 0.], [4, 5, 6, 0., 0.],
                             [7, 8, 9, 0., 0.], [10, 11, 12, 0., 0.]])
        actual = pad_coords_matrices(m1, m2)
        expected = (m1, m2_expected)
        self.assertEqual(actual, expected)

        # second is longer
        m1 = array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        m2 = array(
            [[1, 2, 3, 9], [4, 5, 6, 94], [7, 8, 9, 6], [10, 11, 12, 4]])
        m1_expected = array(
            [[1, 2, 3, 0.], [4, 5, 6, 0.], [7, 8, 9, 0.], [10, 11, 12, 0.]])
        actual = pad_coords_matrices(m1, m2)
        expected = (m1_expected, m2)
        self.assertEqual(actual, expected)

        # equal length, so no change
        m1 = array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        m2 = array([[0, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        actual = pad_coords_matrices(m1, m2)
        expected = (m1, m2)
        self.assertEqual(actual, expected)

    def test_get_procrustes_results(self):
        sample_id_map = {
            'CP3A1': 'S1',
            'CC1A1': 'S2',
            'CC2A1': 'S3',
            'CP1A1': 'S4'}
        actual = get_procrustes_results(self.pcoa1_f, self.pcoa1_f,
                                        sample_id_map=sample_id_map, randomize=None, max_dimensions=None)
        # just some sanity checks as the individual componenets are
        # already tested -- these are based on looking at the output of the
        # run, and testing to ensure that it hasn't changed
        self.assertEqual(
            set(actual[0].split('\n')),
            set('pc vector number\t1\t2\t3\t4\t5\t6\nS1\t0.116737988534\t0.414627960015\t0.201315243115\t0.113769076804\t-0.283025353088\t-0.144278863311\nS2\t-0.238263544222\t-0.37724227779\t-0.169458651217\t0.0305157004776\t0.112181007345\t0.0677415967093\nS3\t-0.199225958574\t-0.250846540029\t-0.119813087305\t-0.155652031006\t0.18495315824\t-0.160875399364\nS4\t0.320751514262\t0.213460857804\t0.0879564954067\t0.0113672537238\t-0.0141088124974\t0.237412665966\n\n\neigvals\t8976580.24393\t6044862.67619\t4372581.39431\t3161360.10319\t2583594.45275\t2407555.39787\n% variation explained\t23.1764657118\t15.6071186064\t11.2894866423\t8.16225689998\t6.67053450426\t6.21602253997'.split('\n')))
        self.assertEqual(
            set(actual[1].split('\n')),
            set('pc vector number\t1\t2\t3\t4\t5\t6\nS1\t0.116737988534\t0.414627960015\t0.201315243115\t0.113769076804\t-0.283025353088\t-0.144278863311\nS2\t-0.238263544222\t-0.37724227779\t-0.169458651217\t0.0305157004776\t0.112181007345\t0.0677415967093\nS3\t-0.199225958574\t-0.250846540029\t-0.119813087305\t-0.155652031006\t0.18495315824\t-0.160875399364\nS4\t0.320751514262\t0.213460857804\t0.0879564954067\t0.0113672537238\t-0.0141088124974\t0.237412665966\n\n\neigvals\t8976580.24393\t6044862.67619\t4372581.39431\t3161360.10319\t2583594.45275\t2407555.39787\n% variation explained\t23.1764657118\t15.6071186064\t11.2894866423\t8.16225689998\t6.67053450426\t6.21602253997'.split('\n')))
        self.assertTrue(actual[2] < 6e-30)

    def test_get_procrustes_results_imprefect_sample_overlap(self):
        sample_id_map = {
            'aaa': 'S0',
            'bbb': 'S1',
            'ccc': 'S2',
            'ddd': 'S3',
            'eee': 'S4'}
        actual = get_procrustes_results(self.pcoa3_f, self.pcoa4_f,
                                        sample_id_map=sample_id_map, randomize=None, max_dimensions=None)
        # Confirm that only the sample ids that are in both procrustes results
        # show up in the output
        for a in actual[:2]:
            self.assertTrue('S1' in a)
            self.assertTrue('S2' in a)
            self.assertTrue('S3' in a)
            self.assertTrue('S0' not in a)
            self.assertTrue('S4' not in a)

    def test_get_procrustes_results_no_sample_overlap(self):
        """ValueError raised on no overlapping sample ids"""
        self.assertRaises(
            ValueError, get_procrustes_results, self.pcoa1_f, self.pcoa3_f,
            sample_id_map=None, randomize=None, max_dimensions=None)

    def test_procrustes_monte_carlo(self):
        """ sanity test of procrustes_monte_carlo wrapper function

          THIS TEST MAY RANDOMLY FAIL BECAUSE IT IS BASED ON TESTING
          RANDOM PERMUTATIONS, BUT THAT SHOULD BE RARE.

        """
        def shuffle_f(coords):
            """ A fake shuffle function -- used to avoid random test failures

                returns a re-ordered coords2
            """
            return array(
                [[-0.16713783, 0.22321481, 0.33766418, 0.22785083, -0.23830446, -0.18754852],
                 [-0.14864343,
                  0.07290181,
                  -0.06250315,
                  0.03201972,
                  -0.0966749,
                  -0.10337987],
                 [0.35725269,
                  -0.00761567,
                  0.09044279,
                  -0.21006839,
                  -0.01355589,
                  -0.04590791],
                    [0.26535811, 0.09772598, 0.04339214, -0.21014987, 0.14089095, -0.10261849]])

        actual = procrustes_monte_carlo(self.pcoa1_f,
                                        self.pcoa2_f,
                                        trials=100,
                                        shuffle_f=shuffle_f)
        # just some sanity checks as the individual componenets are
        # already tested -- these are based on looking at the output of the
        # run, and testing to ensure that it hasn't changed
        expected_actual_m2 = 0.0211
        expected_len_trial_m2 = 100
        expected_count_better = 0
        expected_p_value = 0.0
        self.assertAlmostEqual(actual[0], expected_actual_m2, 3)
        self.assertEqual(len(actual[1]), expected_len_trial_m2)
        self.assertEqual(actual[2], expected_count_better)
        self.assertEqual(actual[3], expected_p_value)


pcoa1_f = """pc vector number	1	2	3	4	5	6
CP3A1	-322.585729836	938.204618621	-28.2137779927	490.569459399	-1046.48732174	-234.500487421
CC1A1	-1190.7802044	-998.399315723	-934.981261143	286.964362952	-79.9676118756	284.018377708
CC2A1	-1095.30958517	-689.284928862	-813.567688811	-168.328870745	98.0045191705	-275.089103439
CP1A1	176.351315302	446.228808693	-305.44480187	240.13477453	-388.822965062	698.967247585


eigvals	8976580.24393	6044862.67619	4372581.39431	3161360.10319	2583594.45275	2407555.39787
% variation explained	23.1764657118	15.6071186064	11.2894866423	8.16225689998	6.67053450426	6.21602253997
"""

pcoa2_f = """pc vector number	1	2	3	4	5	6
CC1A1	0.265358109258	0.0977259770776	0.0433921359913	-0.210149870291	0.140890951431	-0.102618488686
CC2A1	0.357252689812	-0.00761566770389	0.0904427944586	-0.210068386229	-0.0135558941951	-0.0459079084294
CP1A1	-0.148643434206	0.072901814482	-0.0625031523656	0.0320197230867	-0.0966748975994	-0.103379871571
CP3A1	-0.167137833152	0.223214806728	0.337664175642	0.227850833252	-0.238304458085	-0.187548520834


eigvals	1.03654365499	0.486727634877	0.436010533243	0.344489629748	0.325443839964	0.301867574234
% variation explained	17.9408597738	8.42445197874	7.54662266189	5.96254688461	5.63289569997	5.2248294546
"""

pcoa3_f = """pc vector number	1	2	3	4	5	6
aaa	0.265358109258	0.0977259770776	0.0433921359913	-0.210149870291	0.140890951431	-0.102618488686
bbb	0.357252689812	-0.00761566770389	0.0904427944586	-0.210068386229	-0.0135558941951	-0.0459079084294
ccc	-0.148643434206	0.072901814482	-0.0625031523656	0.0320197230867	-0.0966748975994	-0.103379871571
ddd	-0.167137833152	0.223214806728	0.337664175642	0.227850833252	-0.238304458085	-0.187548520834


eigvals	1.03654365499	0.486727634877	0.436010533243	0.344489629748	0.325443839964	0.301867574234
% variation explained	17.9408597738	8.42445197874	7.54662266189	5.96254688461	5.63289569997	5.2248294546
"""

pcoa4_f = """pc vector number	1	2	3	4	5	6
eee	0.265358109258	0.0977259770776	0.0433921359913	-0.210149870291	0.140890951431	-0.102618488686
bbb	0.357252689812	-0.00761566770389	0.0904427944586	-0.210068386229	-0.0135558941951	-0.0459079084294
ccc	-0.148643434206	0.072901814482	-0.0625031523656	0.0320197230867	-0.0966748975994	-0.103379871571
ddd	-0.167137833152	0.223214806728	0.337664175642	0.227850833252	-0.238304458085	-0.187548520834


eigvals	1.03654365499	0.486727634877	0.436010533243	0.344489629748	0.325443839964	0.301867574234
% variation explained	17.9408597738	8.42445197874	7.54662266189	5.96254688461	5.63289569997	5.2248294546
"""

sample_id_map1 = {'CC1A1': 's1', 'CC2A1': 's2', 'CP1A1': 's3', 'CP3A1': 's4',
                  'aaa': 's1', 'bbb': 's2', 'ccc': 's3', 'ddd': 's4'}

if __name__ == "__main__":
    main()
