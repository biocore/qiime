#!/usr/bin/env python

__author__ = "Kyle Patnode"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Kyle Patnode", "Jai Ram Rideout", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Kyle Patnode"
__email__ = "kpatnode1@gmail.com"

"""Test suite for the generate_taxa_compare_table.py module.

Tests each function in the tax2tree controller module. It
should be noted that these tests are fairly sparse, since
tax2tree implements quite a few of its own tests."""

from cogent.app.util import ApplicationNotFoundError
try:
    from t2t.nlevel import load_tree, load_consensus_map, determine_rank_order
except ImportError:
    raise ApplicationNotFoundError(
        "Cannot find tax2tree. Is it installed? Is it in your path?")
from os import makedirs, getcwd, chdir
from os.path import exists
from shutil import rmtree
from tempfile import mkdtemp
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from qiime.test import initiate_timeout, disable_timeout
from qiime.util import get_qiime_temp_dir
from qiime.pycogent_backports.tax2tree import *


class GenerateTaxaCompareTableTests(TestCase):

    """Tests for the tax2tree_controller.py module."""

    def setUp(self):
        """Set up files/environment that will be used by the tests."""
        # The prefix to use for temporary files. This prefix may be added to,
        # but all temp dirs and files created by the tests will have this
        # prefix at a minimum.
        self.prefix = 'tax2tree_controller_tests'

        self.start_dir = getcwd()
        self.dirs_to_remove = []
        self.files_to_remove = []

        self.tmp_dir = get_qiime_temp_dir()
        if not exists(self.tmp_dir):
            makedirs(self.tmp_dir)
            # if test creates the temp dir, also remove it
            self.dirs_to_remove.append(self.tmp_dir)

        initiate_timeout(60)

    def tearDown(self):
        """ """
        disable_timeout()

        # change back to the start dir - some workflows change directory
        chdir(self.start_dir)

        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_clean_output_valid_input(self):
        """Test clean_output using standard valid input."""
        exp = {'e': ('f__Lachnospiraceae; g__Lachnospira; s__', 0),
               'g': ('f__Lachnospiraceae; g__Lachnospira; s__', 0)}
        obs = clean_output(test_results, ['e', 'g'])
        self.assertEquals(obs, exp)

    def test_expand_constrings_valid_input(self):
        """Test expand_constrings using standard valid input

        Makes sure constrings formatting without spaces get them added and
        those with spaces pass through unmolested. Also makes sure 'BAD-
        NAMES' such as Unclassified are handled correctly (e.g. not at all)"""
        exp = dict(x.split('\t') for x in test_cons)
        obs = expand_constrings(test_cons_no_spaces)

        self.assertEqual(obs, exp)

        obs = expand_constrings(test_cons)

        self.assertEqual(obs, exp)

    def test_generate_constrings_valid_input(self):
        """Tests generate_constrings with standard valid input.

        Checks that our output mirrors nlevel (tax2tree's interface)."""
        exp = test_results
        determine_rank_order(test_cons[0].split('\t')[1])
        cons_map = load_consensus_map(test_cons, False)
        tree = load_tree(test_tree, cons_map)

        obs = generate_constrings(tree, cons_map)
        self.assertEqual(obs, exp)

    def test_prep_consensus(self):
        """Tests prep_consensus using standard valid input.

        Simultaneously checks if reference values are being renamed and new
        rep values are being inserted and labeled 'Unclassified'."""
        exp = ['a	f__Lachnospiraceae; g__Bacteroides; s__',
               'b	f__Lachnospiraceae; g__Bacteroides; s__',
               'c	f__Lachnospiraceae; g__Bacteroides; s__Bacteroides pectinophilus',
               'd	Unclassified',
               'd_R	f__Lachnospiraceae; g__Bacteroides; s__Bacteroides pectinophilus',
               'e	Unclassified',
               'f	f__Lachnospiraceae; g__Lachnospira; s__',
               'g	Unclassified',
               'h	Unclassified',
               'h_R	f__Lachnospiraceae; g__Lachnospira; s__Bacteroides pectinophilus',
               'i	Unclassified',
               'j	Unclassified']

        obs = prep_consensus(test_cons, ['d', 'h', 'i', 'j'])

        self.assertEqual(set(obs), set(exp))

# All test data borrowed from Tax2Tree's test code
test_tree = "(((a,b),(c,d)),((e,f),(g,h)));"

test_cons = ['a	f__Lachnospiraceae; g__Bacteroides; s__',
             'b	f__Lachnospiraceae; g__Bacteroides; s__',
             'c	f__Lachnospiraceae; g__Bacteroides; s__Bacteroides pectinophilus',
             'd	f__Lachnospiraceae; g__Bacteroides; s__Bacteroides pectinophilus',
             'e	Unclassified',
             'f	f__Lachnospiraceae; g__Lachnospira; s__',
             'g	Unclassified',
             'h	f__Lachnospiraceae; g__Lachnospira; s__Bacteroides pectinophilus']

test_cons_no_spaces = ['a	f__Lachnospiraceae;g__Bacteroides;s__',
                       'b	f__Lachnospiraceae;g__Bacteroides;s__',
                       'c	f__Lachnospiraceae;g__Bacteroides;s__Bacteroides pectinophilus',
                       'd	f__Lachnospiraceae;g__Bacteroides;s__Bacteroides pectinophilus',
                       'e	Unclassified',
                       'f	f__Lachnospiraceae;g__Lachnospira;s__',
                       'g	Unclassified',
                       'h	f__Lachnospiraceae;g__Lachnospira;s__Bacteroides pectinophilus']

test_results = ['a\tf__Lachnospiraceae; g__Bacteroides; s__',
                'b\tf__Lachnospiraceae; g__Bacteroides; s__',
                'c\tf__Lachnospiraceae; g__Bacteroides; s__Bacteroides pectinophilus',
                'd\tf__Lachnospiraceae; g__Bacteroides; s__Bacteroides pectinophilus',
                'e\tf__Lachnospiraceae; g__Lachnospira; s__',
                'f\tf__Lachnospiraceae; g__Lachnospira; s__',
                'g\tf__Lachnospiraceae; g__Lachnospira; s__',
                'h\tf__Lachnospiraceae; g__Lachnospira; s__']

if __name__ == "__main__":
    main()
