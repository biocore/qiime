from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "1.8.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

from unittest import TestCase, main

from numpy import log, array, nan
from biom import Table

from qiime.index import compute_index


class IndexTests(TestCase):
    def setUp(self):
        self.t1 = Table(array([[0, 1, 2],
                               [0, 0, 1],
                               [1, 1, 0],
                               [3, 0, 1]]),
                        ['O1', 'O2', 'O3', 'O4'],
                        ['S1', 'S2', 'S3'],
                        [{'taxonomy': ['foo', 'bar']},
                         {'taxonomy': ['foo', 'not bar']},
                         {'taxonomy': ['foo', 'bar']},
                         {'taxonomy': ['foo', 'not bar']}])

    def test_compute_index(self):
        increased = set(['bar'])
        decreased = set(['not bar'])

        exp = [('S1', log(1 / 3)), ('S2', nan), ('S3', log(2 / 2))]
        obs = sorted(compute_index(self.t1, increased, decreased))
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
