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
from pandas import DataFrame
import pandas.util.testing as pdt

from skbio.maths.stats.ordination import OrdinationResults

from qiime.trajectory_analysis import normalize_samples


class TrajectoryTests(TestCase):
    """"""
    def setUp(self):
        """"""
        mapdata = {'PC.354': {'BarcodeSequence': 'AGCACGAGCCTA',
                              'Treatment': 'Control',
                              'DOB': '20061218',
                              'Description': '354'},
                   'PC.355': {'BarcodeSequence': 'AACTCGTCGATG',
                              'Treatment': 'Control',
                              'DOB': '20061218',
                              'Description': '355'},
                   'PC.356': {'BarcodeSequence': 'ACAGACCACTCA',
                              'Treatment': 'Control',
                              'DOB': '20061126',
                              'Description': '356'},
                   'PC.481': {'BarcodeSequence': 'ACCAGCGACTAG',
                              'Treatment': 'Control',
                              'DOB': '20070314',
                              'Description': '481'},
                   'PC.593': {'BarcodeSequence': 'AGCAGCACTTGT',
                              'Treatment': 'Control',
                              'DOB': '20071210',
                              'Description': '593'},
                   'PC.607': {'BarcodeSequence': 'AACTGTGCGTAC',
                              'Treatment': 'Fast',
                              'DOB': '20071112',
                              'Description': '607'},
                   'PC.634': {'BarcodeSequence': 'ACAGAGTCGGCT',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '634'},
                   'PC.635': {'BarcodeSequence': 'ACCGCAGAGTCA',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '635'},
                   'PC.636': {'BarcodeSequence': 'ACGGTGAGTGTC',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '636'}}
        self.metamap = DataFrame.from_dict(mapdata, orient='index')

        self.eigvals = np.array([0.512367260461, 0.300719094427,
                                 0.267912066004, 0.208988681078, 0.19169895326,
                                 0.16054234528, 0.15017695712, 0.122457748167,
                                 0.0])
        site = np.array([[-0.258465461183, 0.173999546883, 0.0382875792552,
                          -0.19447750562, 0.0831176020844, 0.262430333201,
                          -0.0231636392235, -0.0184794039581, 0.0],
                         [-0.271001135391, -0.0185951319063, -0.0864841926349,
                          0.118064245315, -0.198808358437, -0.0211723599535,
                          -0.191024027565, 0.155646592377, 0.0],
                         [0.235077898175, 0.0962519254489, -0.345792726714,
                          -0.00320862577619, -0.0963777675519, 0.0457025386953,
                          0.185472813286, 0.0404093971793, 0.0],
                         [0.0261407664325, -0.0111459676533, 0.147660603015,
                          0.29087660853, 0.203945472801, 0.0619712384758,
                          0.101641328709, 0.105690998719, 0.0],
                         [0.285007552283, -0.0192549888483, 0.0623263375385,
                          0.138126799852, -0.104798602423, 0.0951720730628,
                          -0.129636097542, -0.220687170372, 0.0],
                         [0.204636326241, -0.139361150932, 0.291513819623,
                          -0.181566786821, -0.159580132715, -0.0246412130162,
                          0.0866252404441, 0.0996221476871, 0.0],
                         [0.233482403212, 0.225257974068, -0.0188623096268,
                          -0.107729981831, 0.177108999572, -0.192905835151,
                          -0.149819471408, 0.0383549037465, 0.0],
                         [-0.0949631911323, -0.420974802495, -0.154869454869,
                          -0.0898427509281, 0.152618194488, -0.0334232691501,
                          -0.0251224777303, -0.0508988536409, 0.0],
                         [-0.359915158638, 0.113822595435, 0.0662203444138,
                          0.0297579972788, -0.0572254078183, -0.193133506163,
                          0.145026331031, -0.149658611738, 0.0]])
        self.prop_expl = np.array([0.267573832777, 0.15704469605,
                                   0.139911863774, 0.109140272454,
                                   0.100111048503, 0.0838401161912,
                                   0.0784269939011, 0.0639511763509, 0.0])
        site_ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354',
                    'PC.593', 'PC.355', 'PC.607', 'PC.634']
        self.ord_res = OrdinationResults(eigvals=self.eigvals, site=site,
                                         proportion_explained=self.prop_expl,
                                         site_ids=site_ids)

        mapdata = {'PC.355': {'BarcodeSequence': 'AACTCGTCGATG',
                              'Treatment': 'Control',
                              'DOB': '20061218',
                              'Description': '355'},
                   'PC.481': {'BarcodeSequence': 'ACCAGCGACTAG',
                              'Treatment': 'Control',
                              'DOB': '20070314',
                              'Description': '481'},
                   'PC.593': {'BarcodeSequence': 'AGCAGCACTTGT',
                              'Treatment': 'Control',
                              'DOB': '20071210',
                              'Description': '593'},
                   'PC.607': {'BarcodeSequence': 'AACTGTGCGTAC',
                              'Treatment': 'Fast',
                              'DOB': '20071112',
                              'Description': '607'},
                   'PC.635': {'BarcodeSequence': 'ACCGCAGAGTCA',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '635'},
                   'PC.636': {'BarcodeSequence': 'ACGGTGAGTGTC',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '636'}}
        self.subset_metamap = DataFrame.from_dict(mapdata, orient='index')

        site = np.array([[-0.258465461183, 0.173999546883, 0.0382875792552,
                          -0.19447750562, 0.0831176020844, 0.262430333201,
                          -0.0231636392235, -0.0184794039581, 0.0],
                         [-0.271001135391, -0.0185951319063, -0.0864841926349,
                          0.118064245315, -0.198808358437, -0.0211723599535,
                          -0.191024027565, 0.155646592377, 0.0],
                         [0.0261407664325, -0.0111459676533, 0.147660603015,
                          0.29087660853, 0.203945472801, 0.0619712384758,
                          0.101641328709, 0.105690998719, 0.0],
                         [0.204636326241, -0.139361150932, 0.291513819623,
                          -0.181566786821, -0.159580132715, -0.0246412130162,
                          0.0866252404441, 0.0996221476871, 0.0],
                         [0.233482403212, 0.225257974068, -0.0188623096268,
                          -0.107729981831, 0.177108999572, -0.192905835151,
                          -0.149819471408, 0.0383549037465, 0.0],
                         [-0.0949631911323, -0.420974802495, -0.154869454869,
                          -0.0898427509281, 0.152618194488, -0.0334232691501,
                          -0.0251224777303, -0.0508988536409, 0.0]])
        site_ids = ['PC.636', 'PC.635', 'PC.481', 'PC.593', 'PC.355', 'PC.607']
        self.subset_ord_res = OrdinationResults(eigvals=self.eigvals,
                                                site=site,
                                                proportion_explained=
                                                self.prop_expl,
                                                site_ids=site_ids)

        mapdata = {'PC.354': {'BarcodeSequence': 'AGCACGAGCCTA',
                              'Treatment': 'Control',
                              'DOB': '20061218',
                              'Description': '354'},
                   'PC.355': {'BarcodeSequence': 'AACTCGTCGATG',
                              'Treatment': 'Control',
                              'DOB': '20061218',
                              'Description': '355'},
                   'PC.356': {'BarcodeSequence': 'ACAGACCACTCA',
                              'Treatment': 'Control',
                              'DOB': '20061126',
                              'Description': '356'},
                   'PC.481': {'BarcodeSequence': 'ACCAGCGACTAG',
                              'Treatment': 'Control',
                              'DOB': '20070314',
                              'Description': '481'},
                   'PC.607': {'BarcodeSequence': 'AACTGTGCGTAC',
                              'Treatment': 'Fast',
                              'DOB': '20071112',
                              'Description': '607'},
                   'PC.634': {'BarcodeSequence': 'ACAGAGTCGGCT',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '634'},
                   'PC.636': {'BarcodeSequence': 'ACGGTGAGTGTC',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '636'}}
        self.subset_metamap_2 = DataFrame.from_dict(mapdata, orient='index')

        mapdata = {'PC.354': {'BarcodeSequence': 'AGCACGAGCCTA',
                              'Treatment': 'Control',
                              'DOB': '20061218',
                              'Description': '354'},
                   'PC.356': {'BarcodeSequence': 'ACAGACCACTCA',
                              'Treatment': 'Control',
                              'DOB': '20061126',
                              'Description': '356'},
                   'PC.634': {'BarcodeSequence': 'ACAGAGTCGGCT',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '634'}}
        self.metamap_error = DataFrame.from_dict(mapdata, orient='index')

    def test_normalize_samples_equal(self):
        """Does not change any object if the two sets of samples are equal"""
        obs_ord, obs_map = normalize_samples(self.ord_res, self.metamap)

        npt.assert_equal(obs_ord, self.ord_res)
        pdt.assert_frame_equal(obs_map, self.metamap)

    def test_normalize_samples_subset_metamap(self):
        """Returns a subsample of metamap if it has more samples than ord_res
        """
        obs_ord, obs_map = normalize_samples(self.subset_ord_res, self.metamap)

        npt.assert_equal(obs_ord, self.subset_ord_res)
        pdt.assert_frame_equal(obs_map.sort(), self.subset_metamap.sort())

    def test_normalize_samples_subset_ord_res(self):
        """Returns a subsample of ord_res if it has more samples than metamap
        """
        obs_ord, obs_map = normalize_samples(self.ord_res, self.subset_metamap)

        npt.assert_equal(obs_ord, self.subset_ord_res)
        pdt.assert_frame_equal(obs_map, self.subset_metamap)

    def test_normalize_samples_subset_metamap_and_ord_res(self):
        """Returns a subsample of metamap and ord_res with the intersection of
        sample ids"""
        obs_ord, obs_map = normalize_samples(self.subset_ord_res,
                                             self.subset_metamap_2)

        mapdata = {'PC.355': {'BarcodeSequence': 'AACTCGTCGATG',
                              'Treatment': 'Control',
                              'DOB': '20061218',
                              'Description': '355'},
                   'PC.481': {'BarcodeSequence': 'ACCAGCGACTAG',
                              'Treatment': 'Control',
                              'DOB': '20070314',
                              'Description': '481'},
                   'PC.607': {'BarcodeSequence': 'AACTGTGCGTAC',
                              'Treatment': 'Fast',
                              'DOB': '20071112',
                              'Description': '607'},
                   'PC.636': {'BarcodeSequence': 'ACGGTGAGTGTC',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '636'}}
        exp_metamap = DataFrame.from_dict(mapdata, orient='index')

        exp_site = np.array([[-0.258465461183, 0.173999546883, 0.0382875792552,
                              -0.19447750562, 0.0831176020844, 0.262430333201,
                              -0.0231636392235, -0.0184794039581, 0.0],
                             [0.0261407664325, -0.0111459676533,
                              0.147660603015, 0.29087660853, 0.203945472801,
                              0.0619712384758, 0.101641328709, 0.105690998719,
                              0.0],
                             [0.233482403212, 0.225257974068, -0.0188623096268,
                              -0.107729981831, 0.177108999572, -0.192905835151,
                              -0.149819471408, 0.0383549037465, 0.0],
                             [-0.0949631911323, -0.420974802495,
                              -0.154869454869, -0.0898427509281,
                              0.152618194488, -0.0334232691501,
                              -0.0251224777303, -0.0508988536409, 0.0]])
        site_ids = ['PC.636', 'PC.481', 'PC.355', 'PC.607']
        exp_ord_res = OrdinationResults(eigvals=self.eigvals, site=exp_site,
                                        proportion_explained=self.prop_expl,
                                        site_ids=site_ids)

        npt.assert_equal(obs_ord, exp_ord_res)
        pdt.assert_frame_equal(obs_map.sort(), exp_metamap.sort())

    def test_normalize_samples_error(self):
        """Raises an error if `ord_res` and `metamap` has no samples in common
        """
        with self.assertRaises(ValueError):
            normalize_samples(self.subset_ord_res, self.metamap_error)


if __name__ == '__main__':
    main()
