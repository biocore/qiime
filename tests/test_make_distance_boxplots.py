#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Test suite for the make_distance_boxplots.py module."""

from cogent.util.unit_test import TestCase, main
from qiime.make_distance_boxplots import color_field_states

class MakeDistanceBoxplotsTests(TestCase):
    """Tests for the make_distance_boxplots.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.map_f = map_lines.split('\n')
        self.too_many_colors_map_f = too_many_colors_map_lines.split('\n')

    def test_color_field_states(self):
        """Test correctly assigns colors to a field based on another field."""
        # All sample IDs and field states.
        exp = ([(1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)],
                {'y': (0.0, 0.0, 1.0), 'x': (1.0, 0.0, 0.0)})
        obs = color_field_states(self.map_f, ['1', '2', '3', '4', '5', '6'],
                'Foo', ['a', 'b', 'c'], 'Bar')
        self.assertFloatEqual(obs, exp)

        # Subset of sample IDs and field states.
        exp = ([(1.0, 0.0, 0.0)], {'x': (1.0, 0.0, 0.0)})
        obs = color_field_states(self.map_f, ['1', '2'], 'Foo', ['a'], 'Bar')
        self.assertFloatEqual(obs, exp)

        # Color field by itself (useless but still allowed).
        exp = ([(1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.9490196078431372,
                0.45098039215686275, 0.01568627450980392)], {'a':
                (1.0, 0.0, 0.0),
                'c': (0.9490196078431372, 0.45098039215686275,
                0.01568627450980392), 'b': (0.0, 0.0, 1.0)})
        obs = color_field_states(self.map_f, ['1', '2', '3', '4', '5', '6'],
                'Foo', ['a', 'b', 'c'], 'Foo')
        self.assertFloatEqual(obs, exp)

    def test_color_field_states_invalid_input(self):
        """Test correctly raises error on invalid input."""
        # Field to color not in mapping file.
        self.assertRaises(ValueError, color_field_states, self.map_f,
                ['1', '2', '3', '4', '5'], 'Fooz', ['a', 'b'], 'Bar')

        # Field to color by not in mapping file.
        self.assertRaises(ValueError, color_field_states, self.map_f,
                ['1', '2', '3', '4', '5'], 'Foo', ['a', 'b'], 'Barz')

        # Field states are not found in field (due to subset of sample IDs).
        self.assertRaises(ValueError, color_field_states, self.map_f,
                ['1', '2', '3', '4', '5'], 'Foo', ['a', 'c'], 'Bar')

        # Field states are not found in field (not in column at all).
        self.assertRaises(ValueError, color_field_states, self.map_f,
                ['1', '2', '3', '4', '5', '6'], 'Foo', ['a', 'c', 'z'], 'Bar')

        # Not enough colors.
        samp_ids = [str(i) for i in range(1, 31)]
        self.assertRaises(ValueError, color_field_states,
                self.too_many_colors_map_f, samp_ids, 'Description', samp_ids,
                'Description')

        # No one-to-one mapping.
        self.assertRaises(ValueError, color_field_states, self.map_f,
                ['1', '2', '3', '4', '5', '6'], 'Foo', ['a', 'c', 'b'], 'Baz')


map_lines = """#SampleID\tFoo\tBar\tBaz\tDescription
1\ta\tx\tm\t1
2\tb\ty\tn\t2
3\ta\tx\tm\t3
4\ta\tx\tn\t4
5\tb\ty\tn\t5
6\tc\tx\tm\t6"""

too_many_colors_map_lines = """#SampleID\tDescription
1\t1
2\t2
3\t3
4\t4
5\t5
6\t6
7\t7
8\t8
9\t9
10\t10
11\t11
12\t12
13\t13
14\t14
15\t15
16\t16
17\t17
18\t18
19\t19
20\t20
21\t21
22\t22
23\t23
24\t24
25\t25
26\t26
27\t27
28\t28
29\t29
30\t30"""


if __name__ == "__main__":
    main()
