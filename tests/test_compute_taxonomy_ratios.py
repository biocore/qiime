from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "1.9.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

from unittest import TestCase, main

from numpy import log, array, nan, inf
from biom import Table

from qiime.compute_taxonomy_ratios import compute_index


class IndexTests(TestCase):
    def setUp(self):
        ids_and_md = (['O1', 'O2', 'O3', 'O4'],
                      ['S1', 'S2', 'S3'],
                      [{'taxonomy': ['foo', 'bar']},
                       {'taxonomy': ['foo', 'not bar']},
                       {'taxonomy': ['foo', 'bar']},
                       {'taxonomy': ['foo', 'not bar']}])
        self.t1 = Table(array([[0, 1, 2],
                               [0, 0, 1],
                               [1, 1, 0],
                               [3, 0, 1]]), *ids_and_md)
        self.t2 = Table(array([[0, 1, 2],
                               [1, 0, 1],
                               [0, 1, 0],
                               [3, 0, 1]]), *ids_and_md)

    def test_compute_index(self):
        increased = set(['bar'])
        decreased = set(['not bar'])

        exp = [('S1', log(1 / 3)), ('S2', nan), ('S3', log(2 / 2))]
        obs = sorted(compute_index(self.t1, increased, decreased, 'taxonomy'))
        self.assertEqual(obs, exp)

        decreased = set(['foo'])
        increased = set(['bar'])
        exp = [('S1', -inf), ('S2', log(2.0 / 2.0)), ('S3', log(2.0 / 4.0))]
        obs = sorted(compute_index(self.t2, increased, decreased, 'taxonomy'))
        self.assertEqual(obs, exp)

    def test_compute_index_raises(self):
        with self.assertRaises(KeyError):
            next(compute_index(self.t1, set(['foo']), set(['bar']), 'missing'))

        with self.assertRaises(ValueError):
            next(compute_index(self.t1, set(['foo']), set(['x']), 'taxonomy'))


if __name__ == '__main__':
    main()
