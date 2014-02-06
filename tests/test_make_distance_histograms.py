#!/usr/bin/env python
# file test_make_distance_histograms.py

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jeremy Widmann", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"

from cogent.util.unit_test import TestCase, main
import shutil
from os import mkdir, listdir, path
from qiime.parse import parse_mapping_file, parse_distmat, group_by_field,\
    group_by_fields
from numpy import array, arange
from qiime.make_distance_histograms import between_sample_distances, \
    within_category_distances_grouped, between_category_distances_grouped, \
    within_category_distances, within_and_between_fields, \
    all_category_distances, draw_all_histograms, get_histogram_scale, \
    draw_histogram, make_nav_html, make_main_html, get_valid_indices, \
    distances_by_groups, write_distance_files, group_distances, \
    monte_carlo_group_distances, \
    _make_histogram_filenames, _make_relative_paths, \
    _make_random_filename, matplotlib_rgb_color, \
    average_colors, average_all_colors, assign_unassigned_colors,\
    assign_mapped_colors, monte_carlo_group_distances_within_between,\
    get_random_dists, permute_between_groups

from qiime.colors import data_colors
from collections import defaultdict


class DistanceHistogramsTests(TestCase):

    """Tests for make_distance_histograms.py
    """

    def tearDown(self):
        """clean up after running all tests.
        """
        shutil.rmtree(self.working_dir)

    def setUp(self):
        """setup data function for DistanceHistogramsTests."""
        self.working_dir = '/tmp/distance_histogram_tests/'
        try:
            mkdir(self.working_dir)
        except OSError:  # except already exisits
            pass

        self.histogram_dir = path.join(self.working_dir, 'histograms')
        try:
            mkdir(self.histogram_dir)
        except OSError:  # except already exisits remove it and make a new one
            pass

        # Create distance matrix file
        self.dmat_file = self.working_dir + 'dmat.txt'
        dmat_out = open(self.dmat_file, 'w')
        dmat_out.write(DISTANCE_MATRIX_STRING)
        dmat_out.close()

        self.distance_header, self.dmat = \
            parse_distmat(open(self.dmat_file, 'U'))

        # Create mapping file
        self.map_file = self.working_dir + 'map.txt'
        map_out = open(self.map_file, 'w')
        map_out.write(MAPPING_STRING)
        map_out.close()

        mapping, header, comments = parse_mapping_file(
            open(self.map_file, 'U'))
        header[0] = '#' + header[0]
        header = [header]
        header.extend(mapping)
        self.mapping = header

        # Create prefs file
        self.prefs_file = self.working_dir + 'prefs.txt'
        prefs_out = open(self.prefs_file, 'w')
        prefs_out.write(str(PREFS))
        prefs_out.close()

        # Build single field dict for 'Treatment' field.
        self.single_field_treatment = defaultdict(dict)
        self.treatment_groups = group_by_field(self.mapping, 'Treatment')
        self.single_field_treatment['Treatment'] = \
            distances_by_groups(self.distance_header, self.dmat,
                                self.treatment_groups)
        self.paired_field_treatment = {'Treatment_to_Treatment': [
            [('Control', 'Control'), ('Fast', 'Fast'),
             array([[0.729, 0.8, 0.721, 0.765],
                    [0.776, 0.744, 0.749, 0.677],
                    [0.734, 0.777, 0.733, 0.724],
                    [0.696, 0.675, 0.654, 0.696],
                    [0.731, 0.758, 0.738, 0.737]])],
            [('Control', 'Control'), ('Control', 'Control'),
             array([0.625, 0.623, 0.61, 0.577, 0.615,
                    0.642, 0.673, 0.682, 0.737, 0.704])],
            [('Fast', 'Fast'), ('Fast', 'Fast'),
             array([0.718, 0.666, 0.727, 0.6, 0.578, 0.623])]
        ]}

        self.distances_file = self.working_dir + 'distances_out.txt'
        dist_out = open(self.distances_file, 'w')
        dist_out.write(DISTANCES_OUT)
        dist_out.close()

    def test_matplotlib_rgb_color(self):
        """matplotlib_rgb_color should correctly convert RGB to decimal.
        """
        test_data = [(255, 255, 255), (0, 0, 0), (255, 0, 100), (100, 253, 18)]
        exp = [(1.0, 1.0, 1.0), (0., 0., 0.), (1., 0., 0.39215686274509803),
               (0.39215686274509803, 0.99215686274509807, 0.070588235294117646)]
        for t, e in zip(test_data, exp):
            self.assertEqual(matplotlib_rgb_color(t), e)

    def test_average_colors(self):
        """average_colors should properly average two RGB colors.
        """
        to_average = [((255, 255, 255), (0, 0, 0)),
                      ((255, 0, 100), (100, 253, 18))]
        exp = [(127.5, 127.5, 127.5), (177.5, 126.5, 59.)]
        for t, e in zip(to_average, exp):
            self.assertEqual(average_colors(t[0], t[1]), e)

    def test_average_all_colors(self):
        """average_all_colors should properly average all colors.
        """
        to_average = ['Treatment###FIELDDATA###Control_to_Fast',
                      'DOB###FIELDDATA###20070314_to_20061126']
        exp = {'Treatment###FIELDDATA###Control_to_Fast':
               (0.0039215686274509803,
                0.47843137254901963, 0.52745098039215688),
               'DOB###FIELDDATA###20070314_to_20061126': (0.5, 0.0, 0.5)}
        self.assertEqual(average_all_colors(to_average, FIELD_TO_COLOR_PREFS),
                         exp)

    def test_assign_unassigned_colors(self):
        """assign_unassigned_colors should correctly assign unassigned.
        """
        unassigned = ['first', 'second', 'third']
        exp = {'second': (0.6470588235294118, 0.27843137254901962, 0.0),
               'third': (0.9882352941176471, 0.77647058823529413,
                         0.53333333333333333),
               'first': (0.94901960784313721, 0.45098039215686275,
                         0.015686274509803921)}
        self.assertEqual(assign_unassigned_colors(unassigned), exp)

    def test_assign_mapped_colors(self):
        """assign_mapped_colors should correctly assign mapped colors.
        """
        assigned = [
            'Treatment_Within_Control_Distances',
            'DOB_Within_20070314']
        exp = {'Treatment_Within_Control_Distances': (0., 0., 1.),
               'DOB_Within_20070314': (1., 0., 0.)}
        self.assertEqual(assign_mapped_colors(assigned, FIELD_TO_COLOR_PREFS),
                         exp)

    def test_between_sample_distances(self):
        """between_sample_distances should return correct result.
        """
        exp = {'All_Between_Sample_Distances':
               [0.625, 0.623, 0.61, 0.577, 0.729, 0.8, 0.721, 0.765, 0.615,
                0.642, 0.673, 0.776, 0.744, 0.749, 0.677, 0.682, 0.737, 0.734,
                0.777, 0.733, 0.724, 0.704, 0.696, 0.675, 0.654, 0.696, 0.731,
                0.758, 0.738, 0.737, 0.718, 0.666, 0.727, 0.6, 0.578, 0.623]}
        self.assertEqual(between_sample_distances(self.dmat), exp)

    def test_within_category_distances_grouped(self):
        """within_category_distances_grouped should return correct result.
        """
        exp = {'Treatment_All_Within_Category_Distances':
               [0.625, 0.623, 0.61, 0.577, 0.615, 0.642, 0.673, 0.682, 0.737,
                0.704, 0.718, 0.666, 0.727, 0.6, 0.578, 0.623]}
        self.assertEqual(within_category_distances_grouped(
            self.single_field_treatment), exp)

    def test_between_category_distances_grouped(self):
        """between_category_distances_grouped should return correct result.
        """
        exp = {'Treatment_All_Between_Category_Distances':
               [0.729, 0.8, 0.721, 0.765, 0.776, 0.744, 0.749, 0.677,
                0.734, 0.777, 0.733, 0.724, 0.696, 0.675, 0.654, 0.696,
                0.731, 0.758, 0.738, 0.737]}
        self.assertEqual(between_category_distances_grouped(
            self.single_field_treatment), exp)

    def test_within_category_distances(self):
        """within_category_distances should return correct result.
        """
        exp = {'Treatment_Within_Control_Distances':
               [0.625, 0.623, 0.61, 0.577, 0.615,
                0.642, 0.673, 0.682, 0.737, 0.704],
               'Treatment_Within_Fast_Distances':
               [0.718, 0.666, 0.727, 0.6, 0.578, 0.623]}
        self.assertEqual(within_category_distances(
            self.single_field_treatment), exp)

    def test_within_and_between_fields(self):
        """within_and_between_fields should return correct result.
        """
        exp = {'Within_All_Fields':
               [0.625, 0.623, 0.61, 0.577, 0.615, 0.642, 0.673, 0.682, 0.737,
                0.704, 0.718, 0.666, 0.727, 0.6, 0.578, 0.623],
               }

        self.assertEqual(within_and_between_fields(
            self.paired_field_treatment), exp)

    def test_all_category_distances(self):
        """all_category_distances should return correct result.
        """
        exp = {'Treatment###FIELDDATA###Control_to_Control':
               [0.625, 0.623, 0.61, 0.577, 0.615,
                0.642, 0.673, 0.682, 0.737, 0.704],
               'Treatment###FIELDDATA###Control_to_Fast':
               [0.729, 0.8, 0.721, 0.765,
                0.776, 0.744, 0.749, 0.677,
                0.734, 0.777, 0.733, 0.724,
                0.696, 0.675, 0.654, 0.696,
                0.731, 0.758, 0.738, 0.737],
               'Treatment###FIELDDATA###Fast_to_Fast': [0.718, 0.666, 0.727,
                                                        0.6, 0.578, 0.623]
               }
        self.assertEqual(all_category_distances(
            self.single_field_treatment), exp)

    def test_draw_all_histograms(self):
        """draw_all_histograms should return correct result.
        """
        distances_dict, label_to_histogram_filename = \
            draw_all_histograms(single_field=self.single_field_treatment,
                                paired_field=self.paired_field_treatment,
                                dmat=self.dmat,
                                histogram_dir=self.histogram_dir,
                                field_to_color_prefs=FIELD_TO_COLOR_PREFS,
                                background_color='white')

        # Iterate through each histogram file and assure it exisits.
        for k, v in label_to_histogram_filename.items():
            obs_file = open(v, 'U').read()
            self.assertGreaterThan(len(obs_file), 0)

    def test_get_histogram_scale(self):
        """get_histogram_scale should return correct result.
        """
        distances_dict = {'All_Category_Pairs': {'Treatment_Control_to_Control':
                                                 [0.625, 0.623, 0.61, 0.577, 0.615,
                                                  0.642, 0.673, 0.682, 0.737, 0.704],
                                                 'Treatment_Control_to_Fast': [0.729, 0.8, 0.721, 0.765,
                                                                               0.776, 0.744, 0.749, 0.677,
                                                                               0.734, 0.777, 0.733, 0.724,
                                                                               0.696, 0.675, 0.654, 0.696,
                                                                               0.731, 0.758, 0.738, 0.737],
                                                 'Treatment_Fast_to_Fast': [0.718, 0.666, 0.727, 0.6, 0.578,
                                                                            0.623]
                                                 }}
        bins = arange(0, 1.01, 0.05)
        xscale_exp = (0.55, 0.85)
        yscale_exp = (0.0, 0.5)
        xscale_obs, yscale_obs = get_histogram_scale(distances_dict, bins)
        self.assertFloatEqual(xscale_obs, xscale_exp)
        self.assertFloatEqual(yscale_obs, yscale_exp)

    def test_draw_histogram(self):
        """draw_histogram should return correct result.
        """
        hist_outfile = self.working_dir + 'test_hist.png'
        distances = [0.625, 0.623, 0.61, 0.577, 0.615, 0.642, 0.673, 0.682,
                     0.737, 0.704, 0.718, 0.666, 0.727, 0.6, 0.578, 0.623]
        color = 'blue'
        draw_histogram(distances, color, 10, hist_outfile)
        obs_file = open(hist_outfile, 'U').read()
        self.assertGreaterThan(len(obs_file), 0)

    def test_make_nav_html(self):
        """make_nav_html should return correct result.
        """
        distances_dict, label_to_histogram_filename = \
            draw_all_histograms(single_field=self.single_field_treatment,
                                paired_field=self.paired_field_treatment,
                                dmat=self.dmat,
                                histogram_dir=self.histogram_dir,
                                field_to_color_prefs=FIELD_TO_COLOR_PREFS,
                                background_color='white')
        nav_html_obs = \
            make_nav_html(distances_dict, label_to_histogram_filename)
        self.assertEqual(nav_html_obs, NAV_HTML)

    def test_get_valid_indices(self):
        """get_valid_indices should return correct result.
        """
        control_ids = ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593']
        exp_control_indices = [0, 1, 2, 3, 4]

        fast_ids = ['PC.607', 'PC.634', 'PC.635', 'PC.636']
        exp_fast_indices = [5, 6, 7, 8]

        obs_control = get_valid_indices(self.distance_header, control_ids)
        self.assertEqual(obs_control, exp_control_indices)

        obs_fast = get_valid_indices(self.distance_header, fast_ids)
        self.assertEqual(obs_fast, exp_fast_indices)

    def test_distances_by_groups(self):
        """distances_by_groups should return correct result.
        """
        exp = [
            ['Control', 'Fast', array([[0.729, 0.8, 0.721, 0.765],
                                       [0.776, 0.744, 0.749, 0.677],
                                       [0.734, 0.777, 0.733, 0.724],
                                       [0.696, 0.675, 0.654, 0.696],
                                       [0.731, 0.758, 0.738, 0.737]])],
            ['Control', 'Control', array([0.625, 0.623, 0.61, 0.577, 0.615,
                                          0.642, 0.673, 0.682, 0.737, 0.704])],
            ['Fast', 'Fast', array([0.718, 0.666, 0.727, 0.6, 0.578, 0.623])]
        ]

        obs = distances_by_groups(self.distance_header, self.dmat,
                                  self.treatment_groups)
        self.assertEqual(obs, exp)

    def test_write_distance_files(self):
        """write_distance_files should return correct result.
        """
        exp_path = self.working_dir + 'distances/dist_Treatment.txt'
        write_distance_files(self.single_field_treatment,
                             dir_prefix=self.working_dir)
        self.assertEqual(open(exp_path, 'U').read(),
                         open(self.distances_file, 'U').read())

    def test_group_distances(self):
        """group_distances should return correct result.
        """
        single_field_exp = {'Treatment': [
            ['Control', 'Fast', array([[0.729, 0.8, 0.721, 0.765],
                                       [0.776, 0.744, 0.749, 0.677],
                                       [0.734, 0.777, 0.733, 0.724],
                                       [0.696, 0.675, 0.654, 0.696],
                                       [0.731, 0.758, 0.738, 0.737]])],
            ['Control', 'Control', array([0.625, 0.623, 0.61, 0.577, 0.615,
                                          0.642, 0.673, 0.682, 0.737, 0.704])],
            ['Fast', 'Fast', array([0.718, 0.666, 0.727, 0.6, 0.578, 0.623])]
        ]}

        paired_field_exp = {'Treatment_to_Treatment': [
            [('Control', 'Control'), ('Fast', 'Fast'),
             array([[0.729, 0.8, 0.721, 0.765],
                    [0.776, 0.744, 0.749, 0.677],
                    [0.734, 0.777, 0.733, 0.724],
                    [0.696, 0.675, 0.654, 0.696],
                    [0.731, 0.758, 0.738, 0.737]])],
            [('Control', 'Control'), ('Control', 'Control'),
             array([0.625, 0.623, 0.61, 0.577, 0.615,
                    0.642, 0.673, 0.682, 0.737, 0.704])],
            [('Fast', 'Fast'), ('Fast', 'Fast'),
             array([0.718, 0.666, 0.727, 0.6, 0.578, 0.623])]
        ]}

        single_field_obs, paired_field_obs, dmat_obs = \
            group_distances(mapping_file=self.map_file,
                            dmatrix_file=self.dmat_file,
                            fields=['Treatment'],
                            dir_prefix=self.working_dir)

        self.assertEqual(single_field_exp, single_field_obs)
        self.assertEqual(paired_field_exp, paired_field_obs)
        self.assertEqual(self.dmat, dmat_obs)

    def test_monte_carlo_group_distances(self):
        """monte_carlo_group_distances should return correct result.
        """
        mc_group_dist_path = self.working_dir +\
            'monte_carlo_group_distances/group_distances_Treatment.txt'
        monte_carlo_group_distances(mapping_file=self.map_file,
                                    dmatrix_file=self.dmat_file,
                                    prefs=PREFS,
                                    dir_prefix=self.working_dir)
        obs_res = open(mc_group_dist_path, 'U').readlines()
        exp_res = MONTE_CARLO_DISTANCES.split('\n')
        for obs, exp in zip(obs_res, exp_res):
            obs_fields = obs.split('\t')
            exp_fields = exp.split('\t')

            # Check first 8 fields should be identical from run to run.
            for i in range(8):
                self.assertEqual(obs_fields[i], exp_fields[i])

    def test_permute_between_groups(self):
        """permute_between_groups should correctly scrable between groups"""
        a = arange(7)  # 0 - 6
        b = arange(3) + 7  # 7 - 9

        # construct a fake permutation function
        permute_function = lambda size: array([4, 7, 8, 9, 2, 5, 1, 6, 3, 0])
        ar, br = permute_between_groups(a, b, 10, permute_f=permute_function)
        ar_exp = array([4, 7, 8, 9, 2, 5, 1])
        br_exp = array([6, 3, 0])

        # ensure that all permutations are as expected
        for i in xrange(len(ar)):
            self.assertEqual(ar[i], ar_exp)
            self.assertEqual(br[i], br_exp)
        # ensure the correct number of permutations
        self.assertEqual(len(ar), 10)
        self.assertEqual(len(br), 10)

    def test_monte_carlo_group_distances_within_between(self):
        """monte_carlo_group_distances_within_between should return correct.
        """
        mc_group_dist_path = self.working_dir +\
            'monte_carlo_group_distances/group_distances_within_and_between.txt'
        monte_carlo_group_distances_within_between(
            single_field=self.single_field_treatment,
            paired_field=self.paired_field_treatment,
            dmat=self.dmat,
            dir_prefix=self.working_dir)
        obs_res = open(mc_group_dist_path, 'U').readlines()
        exp_res = MONTE_CARLO_DISTANCES_WITHIN_BETWEEN.split('\n')
        for obs, exp in zip(obs_res, exp_res):
            obs_fields = obs.split('\t')
            exp_fields = exp.split('\t')
            # Check first 8 fields should be identical from run to run.
            for i in range(8):
                self.assertEqual(obs_fields[i], exp_fields[i])

    def test_get_random_dists(self):
        """get_random_dists should return correct result.
        """
        real_dists = [['first_a', 'second_a', [.1, .2, .3, .4, .5]],
                      ['first_b', 'second_b', [.01, .02, .03, .04, .05]]]
        obs_dists = get_random_dists(real_dists, self.dmat, 2)
        exp_dists = [[['first_a', 'second_a', [.5, .1, .6, .01, .2]],
                      ['first_b', 'second_b', [.61, .32, .13, .94, .25]],
                      ['first_a', 'second_a', [.7, .02, .12, .40, .05]],
                      ['first_b', 'second_b', [.93, .62, .88, .41, .85]]]]
        for obs, exp in zip(obs_dists[0], exp_dists[0]):
            self.assertEqual(obs[:2], exp[:2])
            self.assertEqual(len(obs[2]), len(exp[2]))

    def test__make_histogram_filenames(self):
        """_make_histogram_filenames should return correct result.
        """
        distances = {'Distances1': [0, .5, .6, .7],
                     'Distances2': [.3, .7, .9, .2, 1.]}
        obs = _make_histogram_filenames(distances, self.histogram_dir)

        self.assertEqual(obs.keys(), distances.keys())
        for k, v in obs.items():
            exp_prefix = self.histogram_dir
            self.assertEquals(path.split(v)[0], exp_prefix)
            self.assertEquals(v[-4:], '.png')

    def test__make_relative_paths(self):
        """_make_relative_paths should return correct result.
        """
        label_to_path = {'first_path': '/path/to/file1.txt',
                         'second_path': '/path/to/file2.txt'}
        prefix = '/path/'
        exp = {'first_path': './to/file1.txt',
               'second_path': './to/file2.txt'}
        self.assertEqual(_make_relative_paths(label_to_path, prefix), exp)

    def test__make_random_filename(self):
        """_make_random_filename should return correct result.
        """
        base_dir = './'
        suffix = 'test_suffix'
        num_chars = 10
        obs = _make_random_filename(base_dir=base_dir, suffix=suffix,
                                    num_chars=num_chars)

        self.assertEqual(
            len(obs), sum([len(base_dir), len(suffix), num_chars]))
        self.assertEqual(obs[:len(base_dir)], base_dir)
        self.assertEqual(obs[len(base_dir) + num_chars:], suffix)


DISTANCE_MATRIX_STRING = \
    """\tPC.354\tPC.355\tPC.356\tPC.481\tPC.593\tPC.607\tPC.634\tPC.635\tPC.636
PC.354\t0.0\t0.625\t0.623\t0.61\t0.577\t0.729\t0.8\t0.721\t0.765
PC.355\t0.625\t0.0\t0.615\t0.642\t0.673\t0.776\t0.744\t0.749\t0.677
PC.356\t0.623\t0.615\t0.0\t0.682\t0.737\t0.734\t0.777\t0.733\t0.724
PC.481\t0.61\t0.642\t0.682\t0.0\t0.704\t0.696\t0.675\t0.654\t0.696
PC.593\t0.577\t0.673\t0.737\t0.704\t0.0\t0.731\t0.758\t0.738\t0.737
PC.607\t0.729\t0.776\t0.734\t0.696\t0.731\t0.0\t0.718\t0.666\t0.727
PC.634\t0.8\t0.744\t0.777\t0.675\t0.758\t0.718\t0.0\t0.6\t0.578
PC.635\t0.721\t0.749\t0.733\t0.654\t0.738\t0.666\t0.6\t0.0\t0.623
PC.636\t0.765\t0.677\t0.724\t0.696\t0.737\t0.727\t0.578\t0.623\t0.0"""

MAPPING_STRING = \
    """#SampleID\tBarcodeSequence\tTreatment\tDOB
PC.354\tAGCACGAGCCTA\tControl\t20061218
PC.355\tAACTCGTCGATG\tControl\t20061218
PC.356\tACAGACCACTCA\tControl\t20061126
PC.481\tACCAGCGACTAG\tControl\t20070314
PC.593\tAGCAGCACTTGT\tControl\t20071210
PC.607\tAACTGTGCGTAC\tFast\t20071112
PC.634\tACAGAGTCGGCT\tFast\t20080116
PC.635\tACCGCAGAGTCA\tFast\t20080116
PC.636\tACGGTGAGTGTC\tFast\t20080116
"""

PREFS = {
    'MONTE_CARLO_GROUP_DISTANCES': {'Treatment': 10},
    'FIELDS': ['Treatment'],
}

DISTANCES_OUT = \
    """Control_to_Fast\t0.729\t0.8\t0.721\t0.765\t0.776\t0.744\t0.749\t0.677\t0.734\t0.777\t0.733\t0.724\t0.696\t0.675\t0.654\t0.696\t0.731\t0.758\t0.738\t0.737
Control_to_Control\t0.625\t0.623\t0.61\t0.577\t0.615\t0.642\t0.673\t0.682\t0.737\t0.704
Fast_to_Fast\t0.718\t0.666\t0.727\t0.6\t0.578\t0.623
"""

MONTE_CARLO_DISTANCES = \
    """Category_1a\tCategory_1b\tAvg\tCategory_2a\tCategory_2b\tAvg\tt\tp\tp_greater\tp_less\tIterations\nControl\tFast\t0.7307\tControl\tControl\t0.6488\t5.13115982149\t1.93767647373e-05\t0.0\t1.0\t10\nControl\tFast\t0.7307\tFast\tFast\t0.652\t3.89776804588\t0.000681957542175\t0.0\t1.0\t10\nControl\tControl\t0.6488\tFast\tFast\t0.652\t-0.11479781274\t0.910235587908\t0.6\t0.4\t10\n"""

MONTE_CARLO_DISTANCES_WITHIN_BETWEEN = \
    """Comparison\tCategory_1\tAvg\tComparison\tCategory_2\tAvg\tt\tp\tp_greater\tp_less\tIterations\nWithin\tTreatment\t0.65\tBetween\tTreatment\t0.7307\t-5.42842519819\t4.76810168062e-06\t1.0\t0.0\t10\nWithin\tTreatment\t0.65\tWithin\tAll_Fields\t0.65\t0.0\t1.0\t0.2\t0.8\t10\nBetween\tTreatment\t0.7307\tWithin\tAll_Fields\t0.65\t5.42842519819\t4.76810168062e-06\t0.0\t1.0\t10\n'"""

NAV_HTML = \
    """<td>\n    <div style="overflow:scroll;white-space:nowrap;width:300px;height: 400px;">\n    <p>\n<span class="normal">All_Within_Category_Grouped</span><br />\n<span class="smnorm"><input type="checkbox" id="check_Treatment_All_Within_Category_Distances"  onclick="visibilityAndOpacity(this, \'Treatment_All_Within_Category_Distances\')" />\n<a onmouseover="mouseoverVisible(\'Treatment_All_Within_Category_Distances\')"; onmouseout="mouseoverHidden(\'Treatment_All_Within_Category_Distances\')">Treatment_All_Within_Category_Distances</a></span><br />\n<span class="normal">All_Between_Category_Grouped</span><br />\n<span class="smnorm"><input type="checkbox" id="check_Treatment_All_Between_Category_Distances"  onclick="visibilityAndOpacity(this, \'Treatment_All_Between_Category_Distances\')" />\n<a onmouseover="mouseoverVisible(\'Treatment_All_Between_Category_Distances\')"; onmouseout="mouseoverHidden(\'Treatment_All_Between_Category_Distances\')">Treatment_All_Between_Category_Distances</a></span><br />\n<span class="normal">All_Between_Sample_Distances</span><br />\n<span class="smnorm"><input type="checkbox" id="check_All_Between_Sample_Distances" checked onclick="visibilityAndOpacity(this, \'All_Between_Sample_Distances\')" />\n<a onmouseover="mouseoverVisible(\'All_Between_Sample_Distances\')"; onmouseout="mouseoverHidden(\'All_Between_Sample_Distances\')">All_Between_Sample_Distances</a></span><br />\n<span class="normal">All_Within_And_Between_Fields</span><br />\n<span class="smnorm"><input type="checkbox" id="check_Within_All_Fields"  onclick="visibilityAndOpacity(this, \'Within_All_Fields\')" />\n<a onmouseover="mouseoverVisible(\'Within_All_Fields\')"; onmouseout="mouseoverHidden(\'Within_All_Fields\')">Within_All_Fields</a></span><br />\n<span class="normal">All_Category_Pairs</span><br />\n<span class="smnorm"><input type="checkbox" id="check_Treatment###FIELDDATA###Control_to_Control"  onclick="visibilityAndOpacity(this, \'Treatment###FIELDDATA###Control_to_Control\')" />\n<a onmouseover="mouseoverVisible(\'Treatment###FIELDDATA###Control_to_Control\')"; onmouseout="mouseoverHidden(\'Treatment###FIELDDATA###Control_to_Control\')">Treatment_Control_to_Control</a></span><br />\n<span class="smnorm"><input type="checkbox" id="check_Treatment###FIELDDATA###Control_to_Fast"  onclick="visibilityAndOpacity(this, \'Treatment###FIELDDATA###Control_to_Fast\')" />\n<a onmouseover="mouseoverVisible(\'Treatment###FIELDDATA###Control_to_Fast\')"; onmouseout="mouseoverHidden(\'Treatment###FIELDDATA###Control_to_Fast\')">Treatment_Control_to_Fast</a></span><br />\n<span class="smnorm"><input type="checkbox" id="check_Treatment###FIELDDATA###Fast_to_Fast"  onclick="visibilityAndOpacity(this, \'Treatment###FIELDDATA###Fast_to_Fast\')" />\n<a onmouseover="mouseoverVisible(\'Treatment###FIELDDATA###Fast_to_Fast\')"; onmouseout="mouseoverHidden(\'Treatment###FIELDDATA###Fast_to_Fast\')">Treatment_Fast_to_Fast</a></span><br />\n<span class="normal">All_Within_Categories</span><br />\n<span class="smnorm"><input type="checkbox" id="check_Treatment_Within_Control_Distances"  onclick="visibilityAndOpacity(this, \'Treatment_Within_Control_Distances\')" />\n<a onmouseover="mouseoverVisible(\'Treatment_Within_Control_Distances\')"; onmouseout="mouseoverHidden(\'Treatment_Within_Control_Distances\')">Treatment_Within_Control_Distances</a></span><br />\n<span class="smnorm"><input type="checkbox" id="check_Treatment_Within_Fast_Distances"  onclick="visibilityAndOpacity(this, \'Treatment_Within_Fast_Distances\')" />\n<a onmouseover="mouseoverVisible(\'Treatment_Within_Fast_Distances\')"; onmouseout="mouseoverHidden(\'Treatment_Within_Fast_Distances\')">Treatment_Within_Fast_Distances</a></span><br />\n    </p>\n    </div>\n    </td>\n\
"""

MAIN_HTML = \
    """
<html><head> <title>
QIIME - Distance Histograms
</title>
<script type="text/javascript" src="./js/histograms.js"></script>
 <style type="text/css">
.smnorm {color: blue; font-family:Arial,Verdana; font-size:10; font-weight: bold;}
.normal {color: black; font-family:Arial,Verdana; font-size:11; font-weight: bold;}

</style>
</head>
<body>
<div id="overDiv" style="position:absolute; visibility:hidden; z-index:1000;"></div>
<div>
    <table>
        <tr>

<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:visible;" id="All_Between_Sample_Distances" name="visible" src="/tmp/distance_histogram_tests/histograms/FSFRIojoGyfL2iGYVFbr.png" border="0"></td>

<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="All_Between_Sample_Distances" name="hidden" src="/tmp/distance_histogram_tests/histograms/FSFRIojoGyfL2iGYVFbr.png" border="0"></td>

<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_Within_Fast_Distances" name="hidden" src="/tmp/distance_histogram_tests/histograms/6aMGfDaBtx1UctPtV4qT.png" border="0"></td>

<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Within_All_Fields" name="hidden" src="/tmp/distance_histogram_tests/histograms/uBnt1Qbg0qoasH8xCg5Q.png" border="0"></td>

<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment###FIELDDATA###Fast_to_Fast" name="hidden" src="/tmp/distance_histogram_tests/histograms/YdQdIIROef1jixUaS3X4.png" border="0"></td>

<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment###FIELDDATA###Control_to_Fast" name="hidden" src="/tmp/distance_histogram_tests/histograms/mPBMbyfSxunyeR2GSQpQ.png" border="0"></td>

<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_All_Between_Category_Distances" name="hidden" src="/tmp/distance_histogram_tests/histograms/3leAJFaXwrMqtJKgTm7v.png" border="0"></td>

<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_All_Within_Category_Distances" name="hidden" src="/tmp/distance_histogram_tests/histograms/JHmUx9dH2lmdmxPlDumP.png" border="0"></td>

<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment###FIELDDATA###Control_to_Control" name="hidden" src="/tmp/distance_histogram_tests/histograms/wfy4pNiN0j7843Tawase.png" border="0"></td>

<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_Within_Control_Distances" name="hidden" src="/tmp/distance_histogram_tests/histograms/geDPwHcu38bQINLYappd.png" border="0"></td>
        </tr>
    </table>
</div>


<table style="position:absolute;left:600">
    <tr>
    <td>
    <div style="overflow:scroll;white-space:nowrap;width:300px;height: 400px;">
    <p>
<span class="normal">All_Within_Category_Grouped</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_All_Within_Category_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_All_Within_Category_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_All_Within_Category_Distances')"; onmouseout="mouseoverHidden('Treatment_All_Within_Category_Distances')">Treatment_All_Within_Category_Distances</a></span><br />
<span class="normal">All_Between_Category_Grouped</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_All_Between_Category_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_All_Between_Category_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_All_Between_Category_Distances')"; onmouseout="mouseoverHidden('Treatment_All_Between_Category_Distances')">Treatment_All_Between_Category_Distances</a></span><br />
<span class="normal">All_Between_Sample_Distances</span><br />
<span class="smnorm"><input type="checkbox" id="check_All_Between_Sample_Distances" checked onclick="visibilityAndOpacity(this, 'All_Between_Sample_Distances')" />
<a onmouseover="mouseoverVisible('All_Between_Sample_Distances')"; onmouseout="mouseoverHidden('All_Between_Sample_Distances')">All_Between_Sample_Distances</a></span><br />
<span class="normal">All_Within_And_Between_Fields</span><br />
<span class="smnorm"><input type="checkbox" id="check_Within_All_Fields"  onclick="visibilityAndOpacity(this, 'Within_All_Fields')" />
<a onmouseover="mouseoverVisible('Within_All_Fields')"; onmouseout="mouseoverHidden('Within_All_Fields')">Within_All_Fields</a></span><br />
<span class="normal">All_Category_Pairs</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment###FIELDDATA###Control_to_Control"  onclick="visibilityAndOpacity(this, 'Treatment###FIELDDATA###Control_to_Control')" />
<a onmouseover="mouseoverVisible('Treatment###FIELDDATA###Control_to_Control')"; onmouseout="mouseoverHidden('Treatment###FIELDDATA###Control_to_Control')">Treatment_Control_to_Control</a></span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment###FIELDDATA###Control_to_Fast"  onclick="visibilityAndOpacity(this, 'Treatment###FIELDDATA###Control_to_Fast')" />
<a onmouseover="mouseoverVisible('Treatment###FIELDDATA###Control_to_Fast')"; onmouseout="mouseoverHidden('Treatment###FIELDDATA###Control_to_Fast')">Treatment_Control_to_Fast</a></span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment###FIELDDATA###Fast_to_Fast"  onclick="visibilityAndOpacity(this, 'Treatment###FIELDDATA###Fast_to_Fast')" />
<a onmouseover="mouseoverVisible('Treatment###FIELDDATA###Fast_to_Fast')"; onmouseout="mouseoverHidden('Treatment###FIELDDATA###Fast_to_Fast')">Treatment_Fast_to_Fast</a></span><br />
<span class="normal">All_Within_Categories</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Within_Control_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_Within_Control_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_Within_Control_Distances')"; onmouseout="mouseoverHidden('Treatment_Within_Control_Distances')">Treatment_Within_Control_Distances</a></span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Within_Fast_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_Within_Fast_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_Within_Fast_Distances')"; onmouseout="mouseoverHidden('Treatment_Within_Fast_Distances')">Treatment_Within_Fast_Distances</a></span><br />
    </p>
    </div>
    </td>

    </tr>
</table>
"""

DATA_COLOR_ORDER = ['blue1', 'lime', 'red1', 'blue2', 'purple1', 'yellow1',
                    'green2', 'red2', 'teal1', 'purple', 'olive', 'silver', 'gray']

FIELD_TO_COLOR_PREFS = {'DOB': (
    {'20070314': ['PC.481'],
     '20071112': ['PC.607'],
     '20080116': ['PC.634',
                  'PC.635',
                  'PC.636'],
     '20061126': ['PC.356'],
     '20061218': ['PC.354',
                  'PC.355'],
     '20071210': ['PC.593']},
    {'20070314': 'red1',
     '20071112': 'blue2',
     '20080116': 'yellow1',
     '20061126': 'blue1',
     '20061218': 'lime',
     '20071210': 'purple1'},
    data_colors,
    DATA_COLOR_ORDER),
    'Treatment': ({'Control': ['PC.354',
                               'PC.355',
                               'PC.356',
                               'PC.481',
                               'PC.593'],
                   'Fast': ['PC.607',
                            'PC.634',
                            'PC.635',
                            'PC.636']},
                  {'Control': 'blue1',
                   'Fast': 'lime'},
                  data_colors,
                  DATA_COLOR_ORDER),
    'BarcodeSequence': ({'ACCAGCGACTAG': ['PC.481'],
                         'ACCGCAGAGTCA': ['PC.635'],
                         'AACTGTGCGTAC': ['PC.607'],
                         'AGCAGCACTTGT': ['PC.593'],
                         'ACAGAGTCGGCT': ['PC.634'],
                         'AACTCGTCGATG': ['PC.355'],
                         'ACGGTGAGTGTC': ['PC.636'],
                         'AGCACGAGCCTA': ['PC.354'],
                         'ACAGACCACTCA': ['PC.356']},
                        {'ACCAGCGACTAG': 'purple1',
                         'ACCGCAGAGTCA': 'yellow1',
                         'AACTGTGCGTAC': 'lime',
                         'ACAGACCACTCA': 'red1',
                         'ACAGAGTCGGCT': 'blue2',
                         'AACTCGTCGATG': 'blue1',
                         'ACGGTGAGTGTC': 'green2',
                         'AGCAGCACTTGT': 'teal1',
                         'AGCACGAGCCTA': 'red2'},
                        data_colors,
                        DATA_COLOR_ORDER),
    'Description': ({'Fasting mouse, I.D. 607': ['PC.607'],
                     'Control mouse, I.D. 481': ['PC.481'],
                     'Control mouse, I.D. 593': ['PC.593'],
                     'Control mouse, I.D. 356': ['PC.356'],
                     'Control mouse, I.D. 354': ['PC.354'],
                     'Control mouse, I.D. 355': ['PC.355'],
                     'Fasting mouse, I.D. 634': ['PC.634'],
                     'Fasting mouse, I.D. 635': ['PC.635'],
                     'Fasting mouse, I.D. 636': ['PC.636']},
                    {'Fasting mouse, I.D. 607': 'yellow1',
                     'Control mouse, I.D. 481': 'blue2',
                     'Control mouse, I.D. 593': 'purple1',
                     'Control mouse, I.D. 356': 'red1',
                     'Control mouse, I.D. 354': 'blue1',
                     'Control mouse, I.D. 355': 'lime',
                     'Fasting mouse, I.D. 634': 'green2',
                     'Fasting mouse, I.D. 635': 'red2',
                     'Fasting mouse, I.D. 636': 'teal1'},
                    data_colors,
                    DATA_COLOR_ORDER),
    'SampleID': ({'PC.636': ['PC.636'],
                  'PC.355': ['PC.355'],
                  'PC.607': ['PC.607'],
                  'PC.634': ['PC.634'],
                  'PC.635': ['PC.635'],
                  'PC.593': ['PC.593'],
                  'PC.356': ['PC.356'],
                  'PC.481': ['PC.481'],
                  'PC.354': ['PC.354']},
                 {'PC.636': 'teal1',
                  'PC.355': 'lime',
                  'PC.607': 'yellow1',
                  'PC.634': 'green2',
                  'PC.635': 'red2',
                  'PC.593': 'purple1',
                  'PC.356': 'red1',
                  'PC.481': 'blue2',
                  'PC.354': 'blue1'},
                 data_colors,
                 DATA_COLOR_ORDER)}

# run tests if called from command line
if __name__ == "__main__":
    main()
