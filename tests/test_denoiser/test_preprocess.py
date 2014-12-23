#!/usr/bin/env python

"""Tests for preprossessing of the denoiser."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from random import sample

from unittest import TestCase, main
from skbio.parse.sequences import parse_fasta
from bfillings.denoiser import FlowgramCollection, Flowgram

from qiime.denoiser.preprocess import sample_mapped_keys, \
    _average_flowgrams, prefix_filter_flowgrams


class Test_preprocess(TestCase):

    def setUp(self):
        self.test_map = {'1': ('a', 'b', 'c'),
                         '2': ('d', 'e', 'f')}

    def test_sample_mapped_keys(self):
        """sample_mapped_keys retuns a dictionary of subsamples."""

        # With num_coverage=1 only the keys will be sampled
        actual = sample_mapped_keys(self.test_map, 1)
        self.assertEqual(actual, {'1': ['1'], '2': ['2']})

        actual = sample_mapped_keys(self.test_map, 3)
        for key in actual.keys():
           # check number of sampled keys
            self.assertEqual(3, len(actual[key]))
            for x in actual[key]:
               # check that sampled key is in the full list
                correct = list(self.test_map[key])
                correct.append(key)
                self.assertTrue(x in correct)

    def test_average_flowgrams(self):
        """_average_flowgrams computes an averaged flowgram for each cluster."""

        fc = FlowgramCollection({'a': '1.0 0.0 0.0 1.0 1.0 1.2 1.2 0.8',
                                 'b': '1.2 1.0 0.0 0.8 1.2 2.4 1.0 0.0'})

        # return the centroid unmodified if sample_mapping = 1
        actual = list(_average_flowgrams({'a': 'b'}, fc, {'a': ['a']}))
        self.assertEqual(actual, [(fc.getFlow('a'), 'a')])

        actual = list(_average_flowgrams({'a': 'b'}, fc, {'a': ['a', 'b']}))
        self.assertEqual(
            actual, [(Flowgram(['1.1 0.5 0.0 0.9 1.1 1.8 1.1 0.4']), 'a')])

    def test_prefix_filter_flowgrams(self):
        """prefix_filter_flowgrams maps all flowgrams which are exact prefixe."""

        fc = FlowgramCollection({'a': '1.0 0.0 0.0 1.0 1.0 1.2 1.2 0.8',
                                 'b': '1.2 1.0 0.0 0.8 1.2 2.4 1.0 0.0',
                                 'c': '2.0 0.0 0.0 1.0',
                                 'd': '1.1 0.3 0.0 1.1'})

        expected_squeeze_map = {'a': ['c', 'd'], 'b': []}
        expected_map = {'a': ['d'], 'b': [], 'c': []}

        obs = prefix_filter_flowgrams(fc)
        self.assertEqual(obs, (3, 4, expected_map))

        obs = prefix_filter_flowgrams(fc, squeeze=True)
        self.assertEqual(obs, (2, 4, expected_squeeze_map))


if __name__ == "__main__":
    main()
