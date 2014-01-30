#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the make_distance_boxplots.py module."""

from cogent.util.unit_test import TestCase, main
from matplotlib.figure import Figure
from qiime.make_distance_boxplots import (_cast_y_axis_extrema,
                                          _color_field_states, make_distance_boxplots,
                                          _sort_distributions_by_median)


class MakeDistanceBoxplotsTests(TestCase):

    """Tests for the make_distance_boxplots.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.map_f = map_lines.split('\n')
        self.dm_f = dm_lines.split('\n')
        self.too_many_colors_map_f = too_many_colors_map_lines.split('\n')

    def test_cast_y_axis_extrema(self):
        """Test correctly assigns colors to a field based on another field."""
        obs = _cast_y_axis_extrema(1.0)
        self.assertFloatEqual(obs, 1.0)

        obs = _cast_y_axis_extrema(1)
        self.assertFloatEqual(obs, 1.0)

        obs = _cast_y_axis_extrema('1.0')
        self.assertFloatEqual(obs, 1.0)

        obs = _cast_y_axis_extrema('1')
        self.assertFloatEqual(obs, 1.0)

        obs = _cast_y_axis_extrema('auto')
        self.assertFloatEqual(obs, None)

    def test_cast_y_axis_extrema_invalid_input(self):
        """Test correctly raises an error on bad input."""
        self.assertRaises(ValueError, _cast_y_axis_extrema, 'foo')

    def test_color_field_states(self):
        """Test correctly assigns colors to a field based on another field."""
        # All sample IDs and field states.
        exp = ([(1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)],
               {'y': (0.0, 0.0, 1.0), 'x': (1.0, 0.0, 0.0)})
        obs = _color_field_states(self.map_f, ['1', '2', '3', '4', '5', '6'],
                                  'Foo', ['a', 'b', 'c'], 'Bar')
        self.assertFloatEqual(obs, exp)

        # Subset of sample IDs and field states.
        exp = ([(1.0, 0.0, 0.0)], {'x': (1.0, 0.0, 0.0)})
        obs = _color_field_states(self.map_f, ['1', '2'], 'Foo', ['a'], 'Bar')
        self.assertFloatEqual(obs, exp)

        # Color field by itself (useless but still allowed).
        exp = ([(1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.9490196078431372,
                0.45098039215686275, 0.01568627450980392)], {'a':
                                                             (1.0, 0.0, 0.0),
                                                             'c': (0.9490196078431372, 0.45098039215686275,
                                                                   0.01568627450980392), 'b': (0.0, 0.0, 1.0)})
        obs = _color_field_states(self.map_f, ['1', '2', '3', '4', '5', '6'],
                                  'Foo', ['a', 'b', 'c'], 'Foo')
        self.assertFloatEqual(obs, exp)

    def test_color_field_states_invalid_input(self):
        """Test correctly raises error on invalid input."""
        # Field to color not in mapping file.
        self.assertRaises(ValueError, _color_field_states, self.map_f,
                          ['1', '2', '3', '4', '5'], 'Fooz', ['a', 'b'], 'Bar')

        # Field to color by not in mapping file.
        self.assertRaises(ValueError, _color_field_states, self.map_f,
                          ['1', '2', '3', '4', '5'], 'Foo', ['a', 'b'], 'Barz')

        # Field states are not found in field (due to subset of sample IDs).
        self.assertRaises(ValueError, _color_field_states, self.map_f,
                          ['1', '2', '3', '4', '5'], 'Foo', ['a', 'c'], 'Bar')

        # Field states are not found in field (not in column at all).
        self.assertRaises(ValueError, _color_field_states, self.map_f,
                          ['1', '2', '3', '4', '5', '6'], 'Foo', ['a', 'c', 'z'], 'Bar')

        # Not enough colors.
        samp_ids = [str(i) for i in range(1, 31)]
        self.assertRaises(ValueError, _color_field_states,
                          self.too_many_colors_map_f, samp_ids, 'Description', samp_ids,
                          'Description')

        # No one-to-one mapping.
        self.assertRaises(ValueError, _color_field_states, self.map_f,
                          ['1', '2', '3', '4', '5', '6'], 'Foo', ['a', 'c', 'b'], 'Baz')

    def test_make_distance_boxplots(self):
        """Test correctly generates plot, raw data, and labels."""
        obs = make_distance_boxplots(self.dm_f, self.map_f, ['Foo', 'Bar'])
        self.assertEqual(len(obs), 2)
        self.assertEqual(obs[0][0], 'Foo')
        self.assertTrue(isinstance(obs[0][1], Figure))
        self.assertEqual(len(obs[0][2]), 7)
        self.assertEqual(len(obs[0][3]), 7)
        self.assertEqual(obs[0][4], [None, None, None, None, None, None, None])

        self.assertEqual(obs[1][0], 'Bar')
        self.assertTrue(isinstance(obs[1][1], Figure))
        self.assertEqual(len(obs[1][2]), 5)
        self.assertEqual(len(obs[1][3]), 5)
        self.assertEqual(obs[1][4], [None, None, None, None, None])

    def test_make_distance_boxplots_suppress_plots(self):
        """Test correctly suppresses different plot types."""
        obs = make_distance_boxplots(self.dm_f, self.map_f, ['Bar'],
                                     suppress_all_within=True)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0][0], 'Bar')
        self.assertTrue(isinstance(obs[0][1], Figure))
        self.assertEqual(len(obs[0][2]), 4)
        self.assertEqual(len(obs[0][3]), 4)
        self.assertEqual(obs[0][4], [None, None, None, None])

        obs = make_distance_boxplots(self.dm_f, self.map_f, ['Bar'],
                                     suppress_all_within=True,
                                     suppress_all_between=True)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0][0], 'Bar')
        self.assertTrue(isinstance(obs[0][1], Figure))
        self.assertEqual(len(obs[0][2]), 3)
        self.assertEqual(len(obs[0][3]), 3)
        self.assertEqual(obs[0][4], [None, None, None])

        obs = make_distance_boxplots(self.dm_f, self.map_f, ['Bar'],
                                     suppress_all_within=True,
                                     suppress_all_between=True,
                                     suppress_individual_within=True)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0][0], 'Bar')
        self.assertTrue(isinstance(obs[0][1], Figure))
        self.assertEqual(len(obs[0][2]), 1)
        self.assertEqual(len(obs[0][3]), 1)
        self.assertEqual(obs[0][4], [None])

        obs = make_distance_boxplots(self.dm_f, self.map_f, ['Bar'],
                                     suppress_all_within=True,
                                     suppress_all_between=True,
                                     suppress_individual_between=True)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0][0], 'Bar')
        self.assertTrue(isinstance(obs[0][1], Figure))
        self.assertEqual(len(obs[0][2]), 2)
        self.assertEqual(len(obs[0][3]), 2)
        self.assertEqual(obs[0][4], [None, None])

        obs = make_distance_boxplots(self.dm_f, self.map_f, ['Bar'],
                                     suppress_individual_within=True,
                                     suppress_individual_between=True)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0][0], 'Bar')
        self.assertTrue(isinstance(obs[0][1], Figure))
        self.assertEqual(len(obs[0][2]), 2)
        self.assertEqual(len(obs[0][3]), 2)
        self.assertEqual(obs[0][4], [None, None])

    def test_make_distance_boxplots_box_color(self):
        """Test correctly colors boxes in a variety of ways."""
        # Single box color for all.
        obs = make_distance_boxplots(self.dm_f, self.map_f, ['Bar'],
                                     box_color='r')
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0][0], 'Bar')
        self.assertTrue(isinstance(obs[0][1], Figure))
        self.assertEqual(len(obs[0][2]), 5)
        self.assertEqual(len(obs[0][3]), 5)
        self.assertEqual(obs[0][4], ['r', 'r', 'r', 'r', 'r'])

        # Single box color, with some plots suppressed.
        obs = make_distance_boxplots(self.dm_f, self.map_f, ['Bar'],
                                     box_color='r',
                                     suppress_individual_within=True)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0][0], 'Bar')
        self.assertTrue(isinstance(obs[0][1], Figure))
        self.assertEqual(len(obs[0][2]), 3)
        self.assertEqual(len(obs[0][3]), 3)
        self.assertEqual(obs[0][4], ['r', 'r', 'r'])

        # Color individual within boxes.
        obs = make_distance_boxplots(self.dm_f, self.map_f, ['Foo'],
                                     color_individual_within_by_field='Bar')
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0][0], 'Foo')
        self.assertTrue(isinstance(obs[0][1], Figure))
        self.assertEqual(len(obs[0][2]), 7)
        self.assertEqual(len(obs[0][3]), 7)
        self.assertEqual(len(obs[0][4]), 7)
        self.assertFloatEqual(obs[0][4], [None, None, (1.0, 0.0, 0.0), (0.0,
                                                                        0.0, 1.0), None, None, None])

        # Color individual within boxes, make sure box_color is ignored.
        obs = make_distance_boxplots(self.dm_f, self.map_f, ['Foo'],
                                     box_color='pink', color_individual_within_by_field='Bar')
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0][0], 'Foo')
        self.assertTrue(isinstance(obs[0][1], Figure))
        self.assertEqual(len(obs[0][2]), 7)
        self.assertEqual(len(obs[0][3]), 7)
        self.assertEqual(len(obs[0][4]), 7)
        self.assertFloatEqual(obs[0][4], [None, None, (1.0, 0.0, 0.0), (0.0,
                                                                        0.0, 1.0), None, None, None])

    def test_make_distance_boxplots_invalid_input(self):
        """Test correctly raises an error on invalid input."""
        # No fields provided.
        self.assertRaises(ValueError, make_distance_boxplots, self.dm_f,
                          self.map_f, None)
        self.assertRaises(ValueError, make_distance_boxplots, self.dm_f,
                          self.map_f, [])

        # Nonexistent field.
        self.assertRaises(ValueError, make_distance_boxplots, self.dm_f,
                          self.map_f, ['Foo', 'foobarbaz'])

        # Invalid width/height.
        self.assertRaises(ValueError, make_distance_boxplots, self.dm_f,
                          self.map_f, ['Foo', 'Bar'], width=-1, height=5)
        self.assertRaises(ValueError, make_distance_boxplots, self.dm_f,
                          self.map_f, ['Foo', 'Bar'], width=1, height=0)

        # Suppress everything.
        self.assertRaises(ValueError, make_distance_boxplots, self.dm_f,
                          self.map_f, ['Foo', 'Bar'], suppress_all_within=True,
                          suppress_all_between=True, suppress_individual_within=True,
                          suppress_individual_between=True)

    def test_sort_distributions_by_median(self):
        """Test correctly sorts distributions by median."""
        exp = ([[0, 0, 0, 1], [2, 1, 1], [1, 2, 3]],
               ['bar', 'baz', 'foo'], ['b', 'r', 'w'])
        obs = _sort_distributions_by_median(
            [[1, 2, 3], [2, 1, 1], [0, 0, 0, 1]], ['foo', 'baz', 'bar'],
            ['w', 'r', 'b'])
        self.assertEqual(obs, exp)


map_lines = """#SampleID\tFoo\tBar\tBaz\tDescription
1\ta\tx\tm\t1
2\tb\ty\tn\t2
3\ta\tx\tm\t3
4\ta\tx\tn\t4
5\tb\ty\tn\t5
6\tc\tx\tm\t6"""

dm_lines = """\t1\t2\t3\t4\t5\t6
1\t0\t1\t2\t4\t7\t11
2\t1\t0\t3\t5\t8\t12
3\t2\t3\t0\t6\t9\t13
4\t4\t5\t6\t0\t10\t14
5\t7\t8\t9\t10\t0\t15
6\t11\t12\t13\t14\t15\t0"""

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
