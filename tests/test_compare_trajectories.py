#!/usr/bin/env python
from __future__ import division

__author__ = "Jose Antonio Navas Molina"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jose Antonio Navas Molina", "Antonio Gonzalez Pena",
               "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jose Antonio Navas Molina"
__email__ = "josenavasmolina@gmail.com"

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd

from skbio.math.stats.ordination import OrdinationResults
from skbio.math.gradient import (GroupResults, CategoryResults,
                                 GradientANOVAResults)

from qiime.compare_trajectories import run_trajectory_analysis


class CompareTrajectoriesTests(TestCase):

    def setUp(self):
        eigvals = np.array([0.512367260461, 0.300719094427, 0.267912066004,
                            0.208988681078, 0.19169895326, 0.16054234528,
                            0.15017695712, 0.122457748167, 0.0])
        site = np.array([[-0.212230626531, 0.216034194368, 0.03532727349,
                          -0.254450494129, -0.0687468542543, 0.231895596562,
                          0.00496549154314, -0.0026246871695,
                          9.73837390723e-10],
                         [-0.277487312135, -0.0295483215975, -0.0744173437992,
                          0.0957182357964, 0.204714844022, -0.0055407341857,
                          -0.190287966833, 0.16307126638, 9.73837390723e-10],
                         [0.220886492631, 0.0874848360559, -0.351990132198,
                          -0.00316535032886, 0.114635191853, -0.00019194106125,
                          0.188557853937, 0.030002427212, 9.73837390723e-10],
                         [0.0308923744062, -0.0446295973489, 0.133996451689,
                          0.29318228566, -0.167812539312, 0.130996149793,
                          0.113551017379, 0.109987942454, 9.73837390723e-10],
                         [0.27616778138, -0.0341866951102, 0.0633000238256,
                          0.100446653327, 0.123802521199, 0.1285839664,
                          -0.132852841046, -0.217514322505, 9.73837390723e-10],
                         [0.202458130052, -0.115216120518, 0.301820871723,
                          -0.18300251046, 0.136208248567, -0.0989435556722,
                          0.0927738484879, 0.0909429797672, 9.73837390723e-10],
                         [0.236467470907, 0.21863434374, -0.0301637746424,
                          -0.0225473129718, -0.205287183891, -0.180224615141,
                          -0.165277751908, 0.0411933458557, 9.73837390723e-10],
                         [-0.105517545144, -0.41405687433, -0.150073017617,
                          -0.116066751485, -0.158763393475, -0.0223918378516,
                          -0.0263068046112, -0.0501209518091,
                          9.73837390723e-10],
                         [-0.371636765565, 0.115484234741, 0.0721996475289,
                          0.0898852445906, 0.0212491652909, -0.184183028843,
                          0.114877153051, -0.164938000185, 9.73837390723e-10]])
        prop_expl = np.array([25.6216900347, 15.7715955926, 14.1215046787,
                              11.6913885817, 9.83044890697, 8.51253468595,
                              7.88775505332, 6.56308246609, 4.42499350906e-16])
        site_ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354', 'PC.593',
                    'PC.355', 'PC.607', 'PC.634']
        self.ord_res = OrdinationResults(eigvals=eigvals, site=site,
                                         proportion_explained=prop_expl,
                                         site_ids=site_ids)
        metadata_map = {'PC.354': {'Treatment': 'Control',
                                   'DOB': '20061218',
                                   'Weight': '60',
                                   'Description': 'Control_mouse_I.D._354'},
                        'PC.355': {'Treatment': 'Control',
                                   'DOB': '20061218',
                                   'Weight': '55',
                                   'Description': 'Control_mouse_I.D._355'},
                        'PC.356': {'Treatment': 'Control',
                                   'DOB': '20061126',
                                   'Weight': '50',
                                   'Description': 'Control_mouse_I.D._356'},
                        'PC.481': {'Treatment': 'Control',
                                   'DOB': '20070314',
                                   'Weight': '52',
                                   'Description': 'Control_mouse_I.D._481'},
                        'PC.593': {'Treatment': 'Control',
                                   'DOB': '20071210',
                                   'Weight': '57',
                                   'Description': 'Control_mouse_I.D._593'},
                        'PC.607': {'Treatment': 'Fast',
                                   'DOB': '20071112',
                                   'Weight': '65',
                                   'Description': 'Fasting_mouse_I.D._607'},
                        'PC.634': {'Treatment': 'Fast',
                                   'DOB': '20080116',
                                   'Weight': '68',
                                   'Description': 'Fasting_mouse_I.D._634'},
                        'PC.635': {'Treatment': 'Fast',
                                   'DOB': '20080116',
                                   'Weight': '70',
                                   'Description': 'Fasting_mouse_I.D._635'},
                        'PC.636': {'Treatment': 'Fast',
                                   'DOB': '20080116',
                                   'Weight': '72',
                                   'Description': 'Fasting_mouse_I.D._636'}}
        self.metadata_map = pd.DataFrame.from_dict(metadata_map,
                                                   orient='index')
        self.categories = ['Treatment']
        self.sort_by = 'Weight'

    # This function makes the comparisons between the results classes easier
    def assert_group_results_almost_equal(self, obs, exp):
        """Tests that obs and exp are almost equal"""
        self.assertEqual(obs.name, exp.name)
        npt.assert_almost_equal(obs.trajectory, exp.trajectory)
        npt.assert_almost_equal(obs.mean, exp.mean)
        self.assertEqual(obs.info.keys(), exp.info.keys())
        for key in obs.info:
            npt.assert_almost_equal(obs.info[key], exp.info[key])
        self.assertEqual(obs.message, exp.message)

    def assert_category_results_almost_equal(self, obs, exp):
        """Tests that obs and exp are almost equal"""
        self.assertEqual(obs.category, exp.category)

        if exp.probability is None:
            self.assertTrue(obs.probability is None)
            self.assertTrue(obs.groups is None)
        else:
            npt.assert_almost_equal(obs.probability, exp.probability)
            for o, e in zip(sorted(obs.groups), sorted(exp.groups)):
                self.assert_group_results_almost_equal(o, e)

    def assert_gradientANOVA_results_almost_equal(self, obs, exp):
        """Tests that obs and exp are almost equal"""
        self.assertEqual(obs.algorithm, exp.algorithm)
        self.assertEqual(obs.weighted, exp.weighted)

        for o, e in zip(sorted(obs.categories), sorted(exp.categories)):
            self.assert_category_results_almost_equal(o, e)

    def test_run_trajectory_analysis_avg(self):
        """Correctly computes the avg method"""
        obs = run_trajectory_analysis(self.ord_res, self.metadata_map,
                                      trajectory_categories=self.categories)
        exp_control_group = GroupResults('Control',
                                         np.array([2.3694943596755276,
                                                   3.3716388181385781,
                                                   5.4452089176253367,
                                                   4.5704258453173559,
                                                   4.4972603724478377]),
                                         4.05080566264,
                                         {'avg': 4.0508056626409275}, None)
        exp_fast_group = GroupResults('Fast', np.array([7.2220488239279126,
                                                        4.2726021564374372,
                                                        1.1169097274372082,
                                                        4.02717600030876]),
                                      4.15968417703,
                                      {'avg': 4.1596841770278292}, None)
        exp_treatment = CategoryResults('Treatment', 0.93311555,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('avg', False, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)

    def test_run_trajectory_analysis_trajectory(self):
        """Correctly computes the trajectory method"""
        obs = run_trajectory_analysis(self.ord_res, self.metadata_map,
                                      trajectory_categories=self.categories,
                                      sort_category=self.sort_by,
                                      algorithm='trajectory')
        exp_control_group = GroupResults('Control', np.array([8.6681963576,
                                                              7.0962717982,
                                                              7.1036434615,
                                                              4.0675712674]),
                                         6.73392072123,
                                         {'trajectory': 13.874494152}, None)
        exp_fast_group = GroupResults('Fast', np.array([11.2291654905,
                                                        3.9163741156,
                                                        4.4943507388]),
                                      6.5466301150,
                                      {'trajectory': 12.713431181}, None)
        exp_treatment = CategoryResults('Treatment', 0.9374500147,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('trajectory', False, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)

    def test_run_trajectory_analysis_diff(self):
        """Correctly computes the first difference method"""
        obs = run_trajectory_analysis(self.ord_res, self.metadata_map,
                                      trajectory_categories=self.categories,
                                      sort_category=self.sort_by,
                                      algorithm='diff')
        exp_control_group = GroupResults('Control', np.array([-1.5719245594,
                                                              0.0073716633,
                                                              -3.0360721941]),
                                         -1.5335416967,
                                         {'mean': -1.5335416967,
                                          'std': 1.2427771485}, None)
        exp_fast_group = GroupResults('Fast', np.array([-7.3127913749,
                                                        0.5779766231]),
                                      -3.3674073758,
                                      {'mean': -3.3674073758,
                                       'std': 3.9453839990}, None)
        exp_treatment = CategoryResults('Treatment', 0.6015260608,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('diff', False, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)

    def test_run_trajectory_analysis_wdiff(self):
        """Correctly computes the window difference method"""
        obs = run_trajectory_analysis(self.ord_res, self.metadata_map,
                                      trajectory_categories=self.categories,
                                      sort_category=self.sort_by,
                                      algorithm='wdiff', window_size=3)
        exp_control_group = GroupResults('Control', np.array([-2.5790341819,
                                                              -2.0166764661,
                                                              -3.0360721941,
                                                              0.]),
                                         -1.9079457105,
                                         {'mean': -1.9079457105,
                                          'std': 1.1592139913}, None)
        exp_fast_group = GroupResults('Fast', np.array([11.2291654905,
                                                        3.9163741156,
                                                        4.4943507388]),
                                      6.5466301150,
                                      {'mean': 6.5466301150,
                                       'std': 3.3194494926},
                                      "Cannot calculate the first difference "
                                      "with a window of size (3).")
        exp_treatment = CategoryResults('Treatment', 0.0103976830,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('wdiff', False, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)

    def test_run_trajectory_analysis_error(self):
        """Raises an error if the algorithm is not recognized"""
        with self.assertRaises(ValueError):
            run_trajectory_analysis(self.ord_res, self.metadata_map,
                                    algorithm='foo')


if __name__ == '__main__':
    main()
