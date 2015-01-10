#!/usr/bin/env python
# file test_plot_taxa_summary.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project"  # consider project name
__credits__ = ["Jesse Stombaugh", "Julia Goodrich"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.0-rc2"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"


import matplotlib
from matplotlib import use
use('Agg', warn=False)
from numpy import array
from os.path import exists
from StringIO import StringIO
from unittest import TestCase, main
from os import remove, mkdir, removedirs, listdir
from qiime.plot_taxa_summary import (make_pie_chart, make_img_name,
                                     get_counts, write_html_file,
                                     make_HTML_table, get_fracs, make_all_charts,
                                     make_area_bar_chart, make_legend)


class TopLevelTests(TestCase):

    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""

        self.props = {"title": "Class: from 3 categories"}
        self.prefs = {'1': {'column': '1'}}
        self.counts1 = [(
            1, "a;b;c", "a<br>b<br>c"), (3, "d;e;f", "d<br>e<br>f"),
            (4, "a;g;h", "a<br>g<br>h"), (2, "d;e;i", "d<br>e<br>i")]
        self.sample_ids = ['14FC041', '14FC042', '14FC043', '14FC044']
        self.taxa = ["a;b;c", "d;e;i", "d;e;f", "a;g;h"]
        self.lines_parsed = (['14FC041', '14FC042', '14FC043', '14FC044'],
                             ['a;b;c', 'd;e;f', 'a;g;h', "d;e;i"],
                             [['0.1', '0.3', '0.2'], ['0', '0.2', '0.1'],
                              ['0.4', '0', '0.3'], ['0.5', '0', '0.1']])
        self.fracs = [("a;b;c", 1.0 / 10), ("d;e;f", 3.0 / 10),
                      ("a;g;h", 4.0 / 10), ("d;e;i", 2.0 / 10)]
        self.colors = ['#0000ff', '#00ff00', '#ff0000', '#00ffff']
        self.area_fracs = [[0.1, 0.3, 0.2], [0.0, 0.2, 0.1],
                           [0.4, 0.0, 0.3], [0.5, 0.0, 0.1]]
        self.color_prefs = {
            "a;b;c": 'blue1', "d;e;i": 'red1', "d;e;f": 'blue2',
            "a;g;h": 'red2'}
        self.dpi = 80
        self.plot_width = 12
        self.plot_height = 6
        self.bar_width = 1
        self.generate_image_type = 'pdf'

        self._paths_to_clean_up = []
        self._dirs_to_clean_up = []
        self.dir_path = "/tmp/qiimewebfiles/"

        # make the webfile directory
        try:
            mkdir(self.dir_path)
        except OSError:
            pass

        # make the charts directory
        try:
            mkdir("/tmp/qiimewebfiles/charts")
        except OSError:
            pass

        # define directory to clean up
        self._dirs_to_clean_up = ["/tmp/qiimewebfiles/charts"]

    def tearDown(self):
        map(remove, self._paths_to_clean_up)
        map(removedirs, self._dirs_to_clean_up)

    def test_make_legend(self):
        """make_legend create a legend image given an array of ids and
            colors"""
        fpath = '/tmp/qiimewebfiles/area.pdf'
        filename1 = '/tmp/qiimewebfiles/area_legend.pdf'

        obs = make_legend(self.sample_ids, self.colors, self.plot_width,
                          self.plot_height, 'black', 'white', fpath,
                          self.generate_image_type, self.dpi)

        self.assertTrue(exists(filename1), 'The png file was not created in \
the appropriate location')

        self._paths_to_clean_up = [filename1]

    def test_get_counts(self):
        """get_counts should gets all the counts for an input file"""

        # test the pie charts
        img_data = get_counts("Phylum", ['14FC041', '14FC042', '14FC043'], 5,
                              "/tmp/qiimewebfiles/", 1, self.lines_parsed,
                              self.prefs, self.color_prefs, 'black', 'white',
                              'pie', self.generate_image_type, self.plot_width,
                              self.plot_height, self.bar_width, self.dpi,
                              'file.txt', 0, 'categorical', False, False)

        self.assertEqual(len(img_data), 8)

        # test the area charts
        img_data = get_counts("Phylum", ['14FC041', '14FC042', '14FC043'], 5,
                              "/tmp/qiimewebfiles/", 1, self.lines_parsed,
                              self.prefs, self.color_prefs, 'black', 'white',
                              'area', self.generate_image_type, self.plot_width,
                              self.plot_height, self.bar_width, self.dpi,
                              'file.txt', 0, 'categorical', False, False)

        self.assertEqual(len(img_data), 2)

        # test the area charts
        img_data = get_counts("Phylum", ['14FC041', '14FC042', '14FC043'], 5,
                              "/tmp/qiimewebfiles/", 1, self.lines_parsed,
                              self.prefs, self.color_prefs, 'black', 'white',
                              'bar', self.generate_image_type, self.plot_width,
                              self.plot_height, self.bar_width, self.dpi,
                              'file.txt', 0, 'categorical', False, False)

        self.assertEqual(len(img_data), 2)

        # clean up files generated
        self._paths_to_clean_up = ["/tmp/qiimewebfiles/charts/" + f
                                   for f in listdir("/tmp/qiimewebfiles/charts")]

    def test_get_fracs(self):
        """"get_fracs should Return fractions for matplotlib chart"""

        # test the pie charts
        exp_all_counts = [DATA_HTML % (
            "4", (4.0 / 10) * 100, 'a<br>g', 'h', 'h', "a;g;h"),
            DATA_HTML % (
                '3',
                (3.0 / 10) * 100,
                'd<br>e',
                'f',
                'f',
                "d;e;f"),
            DATA_HTML % (
                '2',
                (2.0 / 10) * 100,
                'd<br>e',
                'i',
                'i',
                "d;e;i"),
            DATA_HTML % ('1', (1.0 / 10) * 100, 'a<br>b', 'c', 'c', "a;b;c")]

        fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac \
            = get_fracs(self.counts1, 5, 10, 'pie')

        self.assertEqual(
            fracs_labels_other, [("a;b;c", 1.0 / 10), ("a;g;h", 4.0 / 10),
                                 ("d;e;f", 3.0 / 10), ("d;e;i", 2.0 / 10)])
        self.assertEqual(fracs_labels, [])
        self.assertEqual(all_counts, exp_all_counts)

        self.assertEqual(other_cat, 0)
        self.assertEqual(red, 10)
        self.assertEqual(other_frac, 0)

        fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac \
            = get_fracs(self.counts1, 3, 10, 'pie')

        self.assertEqual(
            fracs_labels_other, [("a;g;h", 4.0 / 10), ("d;e;f", 3.0 / 10)])
        self.assertEqual(
            fracs_labels, [("a;g;h", 4.0 / 7), ("d;e;f", 3.0 / 7)])

        self.assertEqual(all_counts, exp_all_counts)
        self.assertEqual(other_cat, 2)
        self.assertEqual(red, 7)
        self.assertEqual(other_frac, 3.0 / 10)

        # test the area charts
        exp_all_counts = ['4', '3', '2', '1']
        fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac \
            = get_fracs(self.counts1, 5, 10, 'area')

        self.assertEqual(
            fracs_labels_other, [("a;b;c", 1.0 / 10), ("a;g;h", 4.0 / 10),
                                 ("d;e;f", 3.0 / 10), ("d;e;i", 2.0 / 10)])
        self.assertEqual(fracs_labels, [])
        self.assertEqual(all_counts, exp_all_counts)

        self.assertEqual(other_cat, 0)
        self.assertEqual(red, 10)
        self.assertEqual(other_frac, 0)

        fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac \
            = get_fracs(self.counts1, 3, 10, 'area')

        self.assertEqual(
            fracs_labels_other, [("a;g;h", 4.0 / 10), ("d;e;f", 3.0 / 10),
                                 ('d;e;i', 2.0 / 10)])
        self.assertEqual(
            fracs_labels, [("a;g;h", 8.0 / 18), ("d;e;f", 6.0 / 18)])

        self.assertEqual(all_counts, exp_all_counts)
        self.assertEqual(other_cat, 2)
        self.assertEqual(red, 9)
        self.assertEqual(other_frac, 3.0 / 10)

        # test bar charts
        exp_all_counts = ['4', '3', '2', '1']
        fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac \
            = get_fracs(self.counts1, 5, 10, 'bar')

        self.assertEqual(
            fracs_labels_other, [("a;b;c", 1.0 / 10), ("a;g;h", 4.0 / 10),
                                 ("d;e;f", 3.0 / 10), ("d;e;i", 2.0 / 10)])
        self.assertEqual(fracs_labels, [])
        self.assertEqual(all_counts, exp_all_counts)

        self.assertEqual(other_cat, 0)
        self.assertEqual(red, 10)
        self.assertEqual(other_frac, 0)

        fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac \
            = get_fracs(self.counts1, 3, 10, 'bar')

        self.assertEqual(
            fracs_labels_other, [("a;g;h", 4.0 / 10), ("d;e;f", 3.0 / 10),
                                 ('d;e;i', 2.0 / 10)])
        self.assertEqual(
            fracs_labels, [("a;g;h", 8.0 / 18), ("d;e;f", 6.0 / 18)])

        self.assertEqual(all_counts, exp_all_counts)
        self.assertEqual(other_cat, 2)
        self.assertEqual(red, 9)
        self.assertEqual(other_frac, 3.0 / 10)

        self._paths_to_clean_up = ["/tmp/qiimewebfiles/charts/" + f
                                   for f in listdir("/tmp/qiimewebfiles/charts")]

    def test_make_HTML_table(self):
        """make_HTML_table should Make HTML tables for one set charts"""

        # test pie charts
        fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac = \
            get_fracs(self.counts1, 5, 10, 'pie')

        img_data = make_HTML_table("Phylum", other_frac, 10, red, other_cat,
                                   fracs_labels_other, fracs_labels,
                                   self.dir_path, all_counts, 1, self.prefs,
                                   self.color_prefs, 'black', 'white', 'pie',
                                   'Test1',
                                   self.generate_image_type, self.plot_width,
                                   self.plot_height, self.bar_width, self.dpi, 0,
                                   'categorical', False, False)

        self.assertEqual(len(img_data), 2)

        # test area charts
        fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac = \
            get_fracs(self.counts1, 5, 10, 'area')

        img_data = make_HTML_table("Phylum", other_frac, 10, red, other_cat,
                                   fracs_labels_other, fracs_labels,
                                   self.dir_path, all_counts, 1, self.prefs,
                                   self.color_prefs, 'black', 'white', 'pie',
                                   'Test1',
                                   self.generate_image_type, self.plot_width,
                                   self.plot_height, self.bar_width, self.dpi, 0,
                                   'categorical', False, False)

        self.assertEqual(len(img_data), 2)

        # test bar charts
        fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac = \
            get_fracs(self.counts1, 5, 10, 'bar')

        img_data = make_HTML_table("Phylum", other_frac, 10, red, other_cat,
                                   fracs_labels_other, fracs_labels,
                                   self.dir_path, all_counts, 1, self.prefs,
                                   self.color_prefs, 'black', 'white', 'pie',
                                   'Test1',
                                   self.generate_image_type, self.plot_width,
                                   self.plot_height, self.bar_width, self.dpi, 0,
                                   'categorical', False, False)

        self.assertEqual(len(img_data), 2)
        self._paths_to_clean_up = ["/tmp/qiimewebfiles/charts/" + f
                                   for f in listdir("/tmp/qiimewebfiles/charts")]

    def test_make_pie_chart(self):
        """make_pie_chart should create HTML source and pdfs for pie_charts"""

        filename1 = '/tmp/qiimewebfiles/charts/pie_chart.png'
        filename2 = '/tmp/qiimewebfiles/charts/pie_chart_legend.pdf'
        filename3 = '/tmp/qiimewebfiles/charts/pie_chart.pdf'

        obs1, obs2, obs3, obs4 = make_pie_chart(self.fracs, self.dir_path, 1,
                                                self.prefs, self.color_prefs, "black", "white",
                                                self.generate_image_type, self.plot_width,
                                                self.plot_height, self.bar_width, self.dpi, False,
                                                file_prefix="pie_chart",
                                                props=self.props)

        self.assertTrue(exists(filename1), 'The png file was not created in \
the appropriate location')
        self.assertTrue(exists(filename2), 'The eps file was not created in \
the appropriate location')
        self.assertTrue(exists(filename3), 'The pdf file was not created in \
the appropriate location')
        self._paths_to_clean_up = ["/tmp/qiimewebfiles/charts/" + f
                                   for f in listdir("/tmp/qiimewebfiles/charts")]

    def test_make_area_bar_chart(self):
        """make_area_bar_chart should create HTML source and pdfs for area
           and bar charts"""

        # following is a list of files being generated
        filename1 = '/tmp/qiimewebfiles/charts/area_chart.png'
        filename2 = '/tmp/qiimewebfiles/charts/area_chart_legend.pdf'
        filename3 = '/tmp/qiimewebfiles/charts/area_chart.pdf'
        filename4 = '/tmp/qiimewebfiles/charts/bar_chart.png'
        filename5 = '/tmp/qiimewebfiles/charts/bar_chart_legend.pdf'
        filename6 = '/tmp/qiimewebfiles/charts/bar_chart.pdf'

        # test area chart
        obs1, obs2, obs3, obs4 = make_area_bar_chart(self.sample_ids,
                                                     self.area_fracs,
                                                     self.taxa, self.dir_path, 1, self.prefs,
                                                     self.color_prefs, "black", "white", "area",
                                                     self.generate_image_type, self.plot_width,
                                                     self.plot_height, self.bar_width,
                                                     self.dpi, 0, 'categorical', False, False,
                                                     "area_chart")

        self.assertTrue(exists(filename1), 'The png file was not created in \
the appropriate location')
        self.assertTrue(exists(filename2), 'The eps file was not created in \
the appropriate location')
        self.assertTrue(exists(filename3), 'The pdf file was not created in \
the appropriate location')

        # test bar chart
        obs1, obs2, obs3, obs4 = make_area_bar_chart(self.sample_ids,
                                                     self.area_fracs,
                                                     self.taxa, self.dir_path, 1, self.prefs,
                                                     self.color_prefs, "black", "white", "bar",
                                                     self.generate_image_type, self.plot_width,
                                                     self.plot_height, self.bar_width, self.dpi,
                                                     0, 'categorical', False, False, "bar_chart",
                                                     self.props)

        self.assertTrue(exists(filename4), 'The png file was not created in \
the appropriate location')
        self.assertTrue(exists(filename5), 'The eps file was not created in \
the appropriate location')
        self.assertTrue(exists(filename6), 'The pdf file was not created in \
the appropriate location')
        self._paths_to_clean_up = ["/tmp/qiimewebfiles/charts/" + f
                                   for f in listdir("/tmp/qiimewebfiles/charts")]

    def test_write_html_file(self):
        "Write html and make sure it gets cleaned up"""

        filename1 = '/tmp/test.html'
        self._paths_to_clean_up = [filename1]

        write_html_file('Test', '/tmp/test.html')

        self.assertTrue(exists(filename1), 'The file was not created in \
the appropriate location')

        self._paths_to_clean_up = [filename1]


# expected results for the unit testing
DATA_HTML = """<tr class=normal><td>%s</td> <td nowrap>%.2f%%</td>\
<td class="normal" ><a onmouseover="return overlib('<b>Taxonomy:</b><br>%s<br>\
<a href=javascript:gg(\\'%s\\');>%s</a> ',STICKY,MOUSEOFF,RIGHT);" \
onmouseout="return nd();">%s</a></td></tr>"""


# run tests if called from command line
if __name__ == "__main__":
    main()
