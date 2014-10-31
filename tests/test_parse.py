#!/usr/bin/env python
# file test_parse.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Justin Kuczynski", "Greg Caporaso",
               "Cathy Lozupone", "Jens Reeder", "Daniel McDonald",
               "Jai Ram Rideout", "Will Van Treuren", "Yoshiki Vazquez-Baeza",
               "Jose Antonio Navas Molina"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os import close
from tempfile import mkstemp

from numpy import array, nan
from StringIO import StringIO
from unittest import TestCase, main
from numpy.testing import assert_almost_equal
from skbio.util import remove_files
from skbio.io import FileFormatError

from qiime.parse import (group_by_field, group_by_fields,
                         parse_distmat, parse_rarefaction_record, parse_rarefaction, parse_coords,
                         parse_classic_otu_table, make_envs_dict, fields_to_dict,
                         parse_rarefaction_fname, parse_qiime_parameters, parse_qiime_config_files,
                         parse_bootstrap_support, parse_distmat_to_dict, parse_taxonomy,
                         parse_mapping_file, parse_metadata_state_descriptions,
                         parse_rarefaction_data, parse_illumina_line, parse_qual_score,
                         parse_qual_scores, QiimeParseError, parse_newick, parse_trflp,
                         parse_taxa_summary_table, parse_prefs_file, parse_mapping_file_to_dict,
                         mapping_file_to_dict, MinimalQualParser, parse_denoiser_mapping,
                         parse_otu_map, parse_sample_id_map, parse_taxonomy_to_otu_metadata,
                         is_casava_v180_or_later, MinimalSamParser)


class TopLevelTests(TestCase):

    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        self.l19_data = array([
            [7, 1, 0, 0, 0, 0, 0, 0, 0],
            [4, 2, 0, 0, 0, 1, 0, 0, 0],
            [2, 4, 0, 0, 0, 1, 0, 0, 0],
            [1, 7, 0, 0, 0, 0, 0, 0, 0],
            [0, 8, 0, 0, 0, 0, 0, 0, 0],
            [0, 7, 1, 0, 0, 0, 0, 0, 0],
            [0, 4, 2, 0, 0, 0, 2, 0, 0],
            [0, 2, 4, 0, 0, 0, 1, 0, 0],
            [0, 1, 7, 0, 0, 0, 0, 0, 0],
            [0, 0, 8, 0, 0, 0, 0, 0, 0],
            [0, 0, 7, 1, 0, 0, 0, 0, 0],
            [0, 0, 4, 2, 0, 0, 0, 3, 0],
            [0, 0, 2, 4, 0, 0, 0, 1, 0],
            [0, 0, 1, 7, 0, 0, 0, 0, 0],
            [0, 0, 0, 8, 0, 0, 0, 0, 0],
            [0, 0, 0, 7, 1, 0, 0, 0, 0],
            [0, 0, 0, 4, 2, 0, 0, 0, 4],
            [0, 0, 0, 2, 4, 0, 0, 0, 1],
            [0, 0, 0, 1, 7, 0, 0, 0, 0]
        ])
        self.l19_sample_names = [
            'sam1', 'sam2', 'sam3', 'sam4', 'sam5', 'sam6',
            'sam7', 'sam8', 'sam9', 'sam_middle', 'sam11', 'sam12', 'sam13',
            'sam14', 'sam15', 'sam16', 'sam17', 'sam18', 'sam19']
        self.l19_taxon_names = ['tax1', 'tax2', 'tax3', 'tax4', 'endbigtaxon',
                                'tax6', 'tax7', 'tax8', 'tax9']

        self.legacy_otu_table1 = legacy_otu_table1
        self.otu_table1 = otu_table1
        self.otu_table_without_leading_comment = \
            otu_table_without_leading_comment
        self.expected_lineages1 = expected_lineages1
        self.taxa_summary1 = taxa_summary1
        self.taxa_summary1_expected = taxa_summary1_expected
        self.otu_table1_floats = otu_table1_floats
        self.files_to_remove = []
        self.denoiser_mapping1 = denoiser_mapping1.split('\n')
        self.sam_data1 = sam_data1.split("\n")
        self.sam1_expected = sam1_expected

    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_MinimalSamParser(self):
        """MinimalSamParser functions as expected"""
        actual = list(MinimalSamParser(self.sam_data1))
        expected = self.sam1_expected
        self.assertEqual(actual, expected)

    def test_is_casava_v180_or_later(self):
        """ is_casava_v180_or_later functions as expected """
        # handles trailing \n
        header_line = "@M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0\n"
        self.assertTrue(is_casava_v180_or_later(header_line))
        # same w no trailing \n
        header_line = "@M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0"
        self.assertTrue(is_casava_v180_or_later(header_line))

        header_line = "@HWUSI-EAS552R_0357:8:1:10040:6364#0/1"
        self.assertFalse(is_casava_v180_or_later(header_line))
        header_line = "@ some misc junk..."
        self.assertFalse(is_casava_v180_or_later(header_line))

        # non-header line raises error
        header_line = "HWUSI-EAS552R_0357:8:1:10040:6364#0/1"
        self.assertRaises(AssertionError, is_casava_v180_or_later, header_line)
        header_line = "M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0"
        self.assertRaises(AssertionError, is_casava_v180_or_later, header_line)

    def test_parse_taxa_summary_table(self):
        """ parse_taxa_summary_table functions as expected """
        actual = parse_taxa_summary_table(self.taxa_summary1.split('\n'))
        self.assertItemsEqual(actual[0], self.taxa_summary1_expected[0])
        self.assertItemsEqual(actual[1], self.taxa_summary1_expected[1])
        assert_almost_equal(actual[2], self.taxa_summary1_expected[2])

    def test_parse_newick(self):
        """parse_newick correctly matches escaped tip names to otu ids
        """
        # confirm that it works without escaped names
        t1 = ('((((tax7:0.1,tax3:0.2):.98,tax8:.3, tax4:.3):.4,'
              '((tax1:0.3, tax6:.09):0.43,tax2:0.4):0.5):.2,'
              '(tax9:0.3, endbigtaxon:.08));')
        expected1 = ['tax7', 'tax3', 'tax8', 'tax4', 'tax1',
                     'tax6', 'tax2', 'tax9', 'endbigtaxon']
        self.assertEqual(set(parse_newick(t1).getTipNames()), set(expected1))
        self.assertEqual(set([tip.Name for tip in parse_newick(t1).tips()]),
                         set(expected1))

        # throw some screwed up names in
        t2 = ('((((tax7:0.1,tax3:0.2):.98,tax8:.3, \'tax4\':.3):.4,'
              "(('ta_______ x1':0.3, tax6:.09):0.43,tax2:0.4):0.5):.2,"
              '(tax9:0.3, endbigtaxon:.08));')
        expected2 = ['tax7', 'tax3', 'tax8', 'tax4', 'ta_______ x1',
                     'tax6', 'tax2', 'tax9', 'endbigtaxon']
        self.assertEqual(set(parse_newick(t2).getTipNames()), set(expected2))
        self.assertEqual(set([tip.Name for tip in parse_newick(t2).tips()]),
                         set(expected2))

    def test_parse_mapping_file(self):
        """parse_mapping_file functions as expected"""
        s1 = ['#sample\ta\tb', '#comment line to skip',
              'x \t y \t z ', ' ', '#more skip', 'i\tj\tk']
        exp = ([['x', 'y', 'z'], ['i', 'j', 'k']],
               ['sample', 'a', 'b'],
               ['comment line to skip', 'more skip'])
        obs = parse_mapping_file(s1)
        self.assertEqual(obs, exp)

        # We don't currently support this, but we should soon...
        # check that first non-comment, non-blank line is used as
        # header
        # s1 = ['sample\ta\tb', '#comment line to skip',\
        # 'x \t y \t z ', ' ', '#more skip', 'i\tj\tk']
        # exp = ([['x','y','z'],['i','j','k']],\
        #        ['sample','a','b'],\
        #        ['comment line to skip','more skip'])
        # obs = parse_mapping_file(s1)
        # self.assertEqual(obs, exp)

        # check that we strip double quotes by default
        s2 = ['#sample\ta\tb', '#comment line to skip',
              '"x "\t" y "\t z ', ' ', '"#more skip"', 'i\t"j"\tk']
        obs = parse_mapping_file(s2)
        self.assertEqual(obs, exp)

    def test_mapping_file_to_dict(self):
        """parse_mapping_file functions as expected"""
        s1 = ['#sample\ta\tb', '#comment line to skip',
              'x \t y \t z ', ' ', '#more skip', 'i\tj\tk']
        exp = ([['x', 'y', 'z'], ['i', 'j', 'k']],
               ['sample', 'a', 'b'],
               ['comment line to skip', 'more skip'])
        mapres = parse_mapping_file(s1)  # map_data, header, comments
        mapdict = mapping_file_to_dict(*mapres[:2])
        expdict = {'x': {'a': 'y', 'b': 'z'}, 'i': {'a': 'j', 'b': 'k'}}
        self.assertEqual(mapdict, expdict)

    def test_parse_mapping_file_to_dict(self):
        """parse_mapping_file functions as expected"""
        s1 = ['#sample\ta\tb', '#comment line to skip',
              'x \t y \t z ', ' ', '#more skip', 'i\tj\tk']
        exp = ([['x', 'y', 'z'], ['i', 'j', 'k']],
               ['sample', 'a', 'b'],
               ['comment line to skip', 'more skip'])
        mapdict, comments = parse_mapping_file_to_dict(s1)
        expdict = {'x': {'a': 'y', 'b': 'z'}, 'i': {'a': 'j', 'b': 'k'}}
        self.assertEqual(mapdict, expdict)
        self.assertEqual(comments, ['comment line to skip', 'more skip'])

    def test_parse_mapping_file_handles_filepath(self):
        """ parse_mapping_file handles being passed a mapping filepath
        """
        fd, fp = mkstemp(prefix='test_parse_mapping_file',
                        suffix='.txt')
        close(fd)
        self.files_to_remove.append(fp)
        open(fp, 'w').write('\n'.join(['#sample\ta\tb',
                                      '#comment line to skip',
                                       'x \t y \t z ', ' ',
                                       '#more skip',
                                       'i\tj\tk']))
        obs = parse_mapping_file(fp)
        exp = ([['x', 'y', 'z'], ['i', 'j', 'k']],
               ['sample', 'a', 'b'],
               ['comment line to skip', 'more skip'])
        self.assertEqual(obs, exp)

    def test_parse_mapping_file_handles_file_handle(self):
        """ parse_mapping_file handles being passed a mapping filepath
        """
        fd, fp = mkstemp(prefix='test_parse_mapping_file',
                        suffix='.txt')
        close(fd)
        self.files_to_remove.append(fp)
        open(fp, 'w').write('\n'.join(['#sample\ta\tb',
                                      '#comment line to skip',
                                       'x \t y \t z ', ' ',
                                       '#more skip',
                                       'i\tj\tk']))
        obs = parse_mapping_file(open(fp))
        exp = ([['x', 'y', 'z'], ['i', 'j', 'k']],
               ['sample', 'a', 'b'],
               ['comment line to skip', 'more skip'])
        self.assertEqual(obs, exp)

    def test_parse_mapping_file_handles_errors(self):
        """parse_mapping_file handles bad mapping files"""
        # Empty file
        self.assertRaises(QiimeParseError,
                          parse_mapping_file,
                          [])
        # string
        self.assertRaises(QiimeParseError,
                          parse_mapping_file,
                          'my_mapping_file.txt')
        # invalid format (no header line with leading # sign)
        self.assertRaises(QiimeParseError,
                          parse_mapping_file,
                          ['sampleID\ta\tb',
                           '1\tf\t43',
                           '2\tt\t44'])
        # invalid format (no non-header lines)
        self.assertRaises(QiimeParseError,
                          parse_mapping_file,
                          ['#sampleID\ta\tb'])
        # invalid format (no header line)
        self.assertRaises(QiimeParseError,
                          parse_mapping_file,
                          ['1\tf\t43',
                           '2\tt\t44'])

    def test_parse_prefs_file(self):
        """parse_prefs_file should correctly eval prefs string.
        """
        # Test good input
        ps1 = """{'bgcolor':'white','colors':
            {'id':'blue','name':'green'},'list':[1,2,3]}"""
        exp1 = {'bgcolor': 'white', 'colors': {'id': 'blue', 'name': 'green'},
                'list': [1, 2, 3]}
        self.assertEqual(parse_prefs_file(ps1), exp1)

        # Test bad input
        # list of valid input rather than multiline string should fail.
        ps_bad_1 = ["{'bgcolor':'white',",
                    "'colors':{'id':'blue','name':'green'}",
                    ",'list':[1,2,3]}"]
        self.assertRaises(QiimeParseError, parse_prefs_file, ps_bad_1)

        # bad data. Can be evaluated but not a dict.
        ps_bad_2 = "[1,2,3]"
        self.assertRaises(QiimeParseError, parse_prefs_file, ps_bad_2)

    def test_group_by_field(self):
        """group_by_field should group table by fields"""
        t = [
            ['#sample', 'loc', 'age'],
            ['a', 'US', '5'],
            ['b', 'US', '10'],
            ['c', 'Mal', '5'],
            ['d', 'Mal', '10'],
            ['e', 'Ven', '5'],
        ]
        self.assertEqual(group_by_field(t, 'loc'),
                         {'US': ['a', 'b'], 'Mal': ['c', 'd'], 'Ven': ['e']})
        self.assertEqual(group_by_field(t, 'age'),
                         {'5': ['a', 'c', 'e'], '10': ['b', 'd']})

    def test_group_by_fields(self):
        """group_by_fields should group table by fields"""
        t = [
            ['#sample', 'loc', 'age', 'mal'],
            ['a', 'US', '5', 'n'],
            ['b', 'US', '10', 'n'],
            ['c', 'Mal', '5', 'y'],
            ['d', 'Mal', '10', 'n'],
            ['e', 'Mal', '5', 'y'],
        ]
        self.assertEqual(group_by_fields(t, ['age', 'loc']),
                         {('5', 'US'): ['a'], ('10', 'US'): ['b'], ('5', 'Mal'): ['c', 'e'],
                          ('10', 'Mal'): ['d']})

    def test_parse_distmat(self):
        """parse_distmat should read distmat correctly"""
        lines = """\ta\tb\tc
a\t0\t1\t2
b\t1\t0\t3.5
c\t1\t3.5\t0
""".splitlines()
        exp = (['a', 'b', 'c'], array([[0, 1, 2], [1, 0, 3.5], [1, 3.5, 0]]))
        obs = parse_distmat(lines)
        self.assertEqual(obs[0], exp[0])
        assert_almost_equal(obs[1], exp[1])

    def test_parse_distmat_to_dict(self):
        """parse_distmat should return dict of distmat"""
        lines = """\ta\tb\tc
a\t0\t1\t2
b\t1\t0\t3.5
c\t1\t3.5\t0
""".splitlines()
        exp = {'a': {'a': 0.0, 'c': 2.0, 'b': 1.0},
               'c': {'a': 1.0, 'c': 0.0, 'b': 3.5},
               'b': {'a': 1.0, 'c': 3.5, 'b': 0.0}}
        obs = parse_distmat_to_dict(lines)
        self.assertEqual(obs, exp)

        # should raise error because row and column headers don't match
        wrong_dist_mat = """\ta\ty\tx
a\t0\t1\t2
b\t1\t0\t3.5
c\t1\t3.5\t0
""".splitlines()
        self.failUnlessRaises(
            AssertionError,
            parse_distmat_to_dict,
            wrong_dist_mat)

    def test_parse_bootstrap_support(self):
        """parse_distmat should read distmat correctly"""
        input_txt = """#\ta\tb\tc.
#more comments here
node2\t0
17node\t0.11922
"""
        lines = input_txt.splitlines()
        exp = {'17node': 0.11922, 'node2': 0.00}
        obs = parse_bootstrap_support(lines)
        self.assertItemsEqual(obs, exp)

    def test_parse_rarefaction_data(self):
        self.data = {}
        self.data['headers'] = ['PD_whole_tree.txt', 'Antibiotics']
        self.data['error'] = {'NA': [0.099969643842700004],
                              'Y':
                              [0.105669693476,
                               1.08546135424,
                               1.5626248357999999],
                              'N': [0.101173002974]}
        self.data['options'] = ['Y', 'NA', 'N']
        self.data['xaxis'] = [10.0, 310.0, 610.0, 910.0, 1210.0, 1510.0,
                              1810.0, 2110.0, 2410.0, 2710.0, 3010.0]
        self.data['series'] = {'NA': [0.88581050485400004],
                               'Y':
                               [0.918845147059,
                                7.1758656176500004,
                                9.9186072941199992],
                               'N': [0.92636763785999998]}
        self.data['color'] = {'NA': '#00ff00', 'Y': '#ff0000', 'N': '#0000ff'}

        self.rarefaction_series_data = ['# PD_whole_tree.txt',
                                        '# Antibiotics',
                                        'xaxis: 10.0\t310.0\t610.0\t910.0\t1210.0\t1510.0\t1810.0\t2110.0\
        \t2410.0\t2710.0\t3010.0\t',
                                        'xmax: 3310.0',
                                        '>> Y',
                                        'color #ff0000',
                                        'series 0.918845147059\t7.17586561765\t9.91860729412\t',
                                        'error 0.105669693476\t1.08546135424\t1.5626248358\t',
                                        '>> NA',
                                        'color #00ff00',
                                        'series 0.885810504854\t',
                                        'error 0.0999696438427\t',
                                        '>> N',
                                        'color #0000ff',
                                        'series 0.92636763786\t',
                                        'error 0.101173002974'
                                        ]
        test = parse_rarefaction_data(self.rarefaction_series_data)
        self.assertEqual(test, self.data)

    def test_parse_rarefaction_record(self):
        self.rarefactionline1 = 'rare10.txt\t10\t0\t1.99181\t0.42877\t2.13996'
        test1 = parse_rarefaction_record(self.rarefactionline1)
        self.rarefactiondata1 = ('rare10.txt', [10.0, 0.0,
                                                1.9918100000000001, 0.42876999999999998, 2.1399599999999999])
        self.assertEqual(self.rarefactiondata1, test1)

        self.rarefactionline2 = 'rare10.txt\t10\t0\t1.99181\t0.42877\tNA'
        test2 = parse_rarefaction_record(self.rarefactionline2)
        self.rarefactiondata2 = ('rare10.txt', [10.0, 0.0, 1.9918100000000001,
                                                0.42876999999999998, nan])
        self.assertEqual(self.rarefactiondata2, test2)

    def test_parse_rarefaction_fname(self):
        """ parse_rarefaction_fname should return base, seqs/sam, iters, etc."""
        fname = "alpha_rarefaction_900_3.txt"
        base, seqs, iter, ext = parse_rarefaction_fname(fname)
        self.assertEqual((base, seqs, iter, ext),
                         ("alpha_rarefaction", 900, 3, ".txt"))

    def test_parse_rarefaction(self):
        self.rarefactionfile = [
            '\tsequences per sample\titeration\t123\t234\t345',
            'rare10.txt\t10\t0\t1.99181\t0.42877\t2.13996',
            'rare10.txt\t10\t1\t2.07163\t0.42877\t2.37055',
            'rare310.txt\t310\t0\t8.83115\t0.42877\t11.00725',
            'rare310.txt\t310\t1\t10.05242\t0.42877\t8.24474',
            'rare610.txt\t610\t0\t12.03067\t0.42877\t11.58928',
            'rare610.txt\t610\t1\t12.9862\t0.42877\t11.58642']

        self.col_headers = [
            '',
            'sequences per sample',
            'iteration',
            '123',
            '234',
            '345']
        self.comments = []
        self.rarefaction_fns = [
            'rare10.txt',
            'rare10.txt',
            'rare310.txt',
            'rare310.txt',
            'rare610.txt',
            'rare610.txt']
        self.rarefaction_data = [[10.0,
                                  0.0,
                                  1.9918100000000001,
                                  0.42876999999999998,
                                  2.1399599999999999],
                                 [10.0,
                                  1.0,
                                  2.0716299999999999,
                                  0.42876999999999998,
                                  2.3705500000000002],
                                 [310.0,
                                  0.0,
                                  8.8311499999999992,
                                  0.42876999999999998,
                                  11.007250000000001],
                                 [310.0,
                                  1.0,
                                  10.05242,
                                  0.42876999999999998,
                                  8.2447400000000002],
                                 [610.0,
                                  0.0,
                                  12.030670000000001,
                                  0.42876999999999998,
                                  11.58928],
                                 [610.0,
                                  1.0,
                                  12.9862,
                                  0.42876999999999998,
                                  11.58642]]

        test_col_headers, test_comments, test_rarefaction_fns, test_rarefaction_data = parse_rarefaction(
            self.rarefactionfile)
        self.assertEqual(test_col_headers, self.col_headers)
        self.assertEqual(test_comments, self.comments)
        self.assertEqual(test_rarefaction_fns, self.rarefaction_fns)
        self.assertEqual(test_rarefaction_data, self.rarefaction_data)

#     def test_parse_rarefaction(self):
#         """parse_rarefaction should handle multiple recs"""
# recs ="""#HEADER  97.0    NFkeyRightShift 1000
# CHAO1    288.12903   241.27093   371.41733
# ACE  294.74813   252.75930   361.30606   0.68654
# SHANNON  3.71822 3.61146 3.82499
# SIMPSON  0.08021
# n    rare    rare_lci    rare_hci
# 1 1.000000    1.000000    1.000000
# 51    26.878000   26.032290   27.723710
# HEADER   97.0    DMkeySpace  1000
# CHAO1    90.20000    53.81120    196.98011
# ACE  122.66234   68.48901    264.46888   1.53571
# SHANNON  1.35156 1.12549 1.57762
# SIMPSON  0.56002
# n    rare    rare_lci    rare_hci
# 1 1.000000    1.000000    1.000000
# 51    10.707000   10.085986   11.328014
# 101   17.547000   17.046632   18.047368
# 151   23.410000   23.009805   23.810195
# HEADER   97.0    RKkeyW  1000
# CHAO1    251.57143   176.86438   401.83780
# ACE  264.44294   192.23364   395.08515   0.89714
# SHANNON  1.56521 1.44411 1.68630
# SIMPSON  0.54772
# n    rare    rare_lci    rare_hci
# 1 1.000000    1.000000    1.000000
# 51    11.739000   11.060658   12.417342
# 101   19.054000   18.458137   19.649863
# 151   24.994000   24.447174   25.540826
# """.splitlines()
#         obs = parse_rarefaction(recs)
#         self.assertEqual(set(obs.keys()), \
#             set(['NFkeyRightShift','DMkeySpace','RKkeyW']))

    def test_parse_coords(self):
        """parse_coords should handle coords file"""
        coords = ["Eigvals\t3",
                  "4.94\t1.79\t1.50",
                  "",
                  "Proportion explained\t3",
                  "14.3\t5.2\t4.3",
                  "",
                  "Species\t0\t0",
                  "",
                  "Site\t3\t3",
                  "A\t.11\t.09\t.23",
                  "B\t.03\t.07\t-.26",
                  "C\t.12\t.06\t-.32",
                  "",
                  "Biplot\t0\t0",
                  "",
                  "Site constraints\t0\t0"]

        obs = parse_coords(coords)

        exp = (['A', 'B', 'C'],
               array([[.11, .09, .23], [.03, .07, -.26], [.12, .06, -.32]]),
               array([4.94, 1.79, 1.50]),
               array([14.3, 5.2, 4.3]))
        self.assertEqual(obs[0], exp[0])
        assert_almost_equal(obs[1], exp[1])

    def test_parse_coords_exceptions(self):
        """Check exceptions are raised accordingly with missing information"""

        # missing eigenvalues line
        with self.assertRaises(FileFormatError):
            out = parse_coords(COORDS_NO_EIGENVALS.splitlines())
        # missing percentages explained line
        with self.assertRaises(FileFormatError):
            out = parse_coords(COORDS_NO_PCNTS.splitlines())
        # missing vector number line
        with self.assertRaises(FileFormatError):
            out = parse_coords(COORDS_NO_VECTORS.splitlines())

        # a whole different file (taxa summary)
        with self.assertRaises(FileFormatError):
            out = parse_coords(taxa_summary1.splitlines())

    def test_parse_classic_otu_table_legacy(self):
        """parse_classic_otu_table functions as expected with legacy OTU table
        """
        data = self.legacy_otu_table1
        data_f = (data.split('\n'))
        obs = parse_classic_otu_table(data_f)
        exp = (['Fing', 'Key', 'NA'],
               ['0', '1', '2', '3', '4'],
               array([[19111, 44536, 42], [1216, 3500, 6], [1803, 1184, 2],
                      [1722, 4903, 17], [589, 2074, 34]]),
               self.expected_lineages1)
        self.assertEqual(obs[0], exp[0])
        self.assertEqual(obs[1], exp[1])
        assert_almost_equal(obs[2], exp[2])
        self.assertEqual(obs[3], exp[3])

    def test_parse_classic_otu_table(self):
        """parse_classic_otu_table functions as expected with new-style OTU table
        """
        data = self.otu_table1
        data_f = (data.split('\n'))
        obs = parse_classic_otu_table(data_f)
        exp = (['Fing', 'Key', 'NA'],
               ['0', '1', '2', '3', '4'],
               array([[19111, 44536, 42], [1216, 3500, 6], [1803, 1184, 2],
                      [1722, 4903, 17], [589, 2074, 34]]),
               self.expected_lineages1)
        self.assertEqual(obs[0], exp[0])
        self.assertEqual(obs[1], exp[1])
        assert_almost_equal(obs[2], exp[2])
        self.assertEqual(obs[3], exp[3])


        # test that the modified parse_classic performs correctly on OTU tables
        # without leading comments
        data = self.otu_table_without_leading_comment
        data_f = (data.split('\n'))
        obs = parse_classic_otu_table(data_f)
        sams = ['let-7i', 'miR-7', 'miR-17n', 'miR-18a', 'miR-19a', 'miR-22',
                'miR-25', 'miR-26a']
        otus = ['A2M', 'AAAS', 'AACS', 'AADACL1']
        vals = array([
            [-0.2, 0.03680505, 0.205, 0.23, 0.66, 0.08, -0.373, 0.26],
            [-0.09, -0.25, 0.274, 0.15, 0.12, 0.29, 0.029, -0.1148452],
            [0.33, 0.19, 0.27, 0.28, 0.19, 0.25, 0.089, 0.14],
            [0.49, -0.92, -0.723, -0.23, 0.08, 0.49, -0.386, -0.64]])
        exp = (sams, otus, vals, [])  # no lineages
        # because float comps in arrays always errors
        self.assertEqual(obs[0], exp[0])
        self.assertEqual(obs[1], exp[1])
        self.assertEqual(obs[3], exp[3])
        self.assertTrue(all((obs[2] == exp[2]).tolist()))

    def test_parse_classic_otu_table_floats_in_table(self):
        """parse_classic_otu_table functions using an OTU table containing floats
           but cast as int....this will automatically cast into floats"""

        data = self.otu_table1_floats
        data_f = (data.split('\n'))
        obs = parse_classic_otu_table(data_f)
        exp = (['Fing', 'Key', 'NA'],
               ['0', '1', '2', '3', '4'],
               array([[19111.0, 44536.0, 42.0], [1216.0, 3500.0, 6.0],
                      [1803.0, 1184.0, 2.0], [1722.1, 4903.2, 17.0],
                      [589.6, 2074.4, 34.5]]),
               self.expected_lineages1)
        self.assertEqual(obs[0], exp[0])
        self.assertEqual(obs[1], exp[1])
        assert_almost_equal(obs[2], exp[2])
        self.assertEqual(obs[3], exp[3])


    def test_parse_classic_otu_table_float_counts(self):
        """parse_classic_otu_table should return correct result from small table"""
        data = """#Full OTU Counts
#OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors"""
        data_f = (data.split('\n'))
        obs = parse_classic_otu_table(data_f, count_map_f=float)
        exp = (['Fing', 'Key', 'NA'],
               ['0', '1', '2', '3', '4'],
               array(
                   [[19111., 44536., 42.], [1216., 3500., 6.], [1803., 1184., 2.],
                    [1722., 4903., 17.], [589, 2074., 34.]]),
               [['Bacteria', 'Actinobacteria', 'Actinobacteridae', 'Propionibacterineae', 'Propionibacterium'],
                ['Bacteria',
                 'Firmicutes',
                 'Alicyclobacillaceae',
                 'Bacilli',
                 'Lactobacillales',
                 'Lactobacillales',
                 'Streptococcaceae',
                 'Streptococcus'],
                ['Bacteria',
                 'Actinobacteria',
                 'Actinobacteridae',
                 'Gordoniaceae',
                 'Corynebacteriaceae'],
                ['Bacteria',
                 'Firmicutes',
                 'Alicyclobacillaceae',
                 'Bacilli',
                 'Staphylococcaceae'],
                ['Bacteria', 'Cyanobacteria', 'Chloroplasts', 'vectors']])
        self.assertEqual(obs[0], exp[0])
        self.assertEqual(obs[1], exp[1])
        assert_almost_equal(obs[2], exp[2])
        self.assertEqual(obs[3], exp[3])


    def test_parse_classic_otu_table_file(self):
        """parse_classic_otu_table should return correct result on fileio format object"""
        data = """#Full OTU Counts
#OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors"""
        data_f = StringIO(data)
        obs = parse_classic_otu_table(data_f)
        exp = (['Fing', 'Key', 'NA'],
               ['0', '1', '2', '3', '4'],
               array([[19111, 44536, 42], [1216, 3500, 6], [1803, 1184, 2],
                      [1722, 4903, 17], [589, 2074, 34]]),
               [['Bacteria', 'Actinobacteria', 'Actinobacteridae', 'Propionibacterineae', 'Propionibacterium'],
                ['Bacteria',
                 'Firmicutes',
                 'Alicyclobacillaceae',
                 'Bacilli',
                 'Lactobacillales',
                 'Lactobacillales',
                 'Streptococcaceae',
                 'Streptococcus'],
                ['Bacteria',
                 'Actinobacteria',
                 'Actinobacteridae',
                 'Gordoniaceae',
                 'Corynebacteriaceae'],
                ['Bacteria',
                 'Firmicutes',
                 'Alicyclobacillaceae',
                 'Bacilli',
                 'Staphylococcaceae'],
                ['Bacteria', 'Cyanobacteria', 'Chloroplasts', 'vectors']])
        self.assertEqual(obs[0], exp[0])
        self.assertEqual(obs[1], exp[1])
        assert_almost_equal(obs[2], exp[2])
        self.assertEqual(obs[3], exp[3])


    def test_parse_classic_otu_table_consensus_lineage(self):
        """parse_classic_otu_table should accept 'consensusLineage'"""
        data = """#Full OTU Counts
#OTU ID	Fing	Key	NA	consensusLineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors"""
        data_f = StringIO(data)
        obs = parse_classic_otu_table(data_f)
        exp = (['Fing', 'Key', 'NA'],
               ['0', '1', '2', '3', '4'],
               array([[19111, 44536, 42], [1216, 3500, 6], [1803, 1184, 2],
                      [1722, 4903, 17], [589, 2074, 34]]),
               [['Bacteria', 'Actinobacteria', 'Actinobacteridae', 'Propionibacterineae', 'Propionibacterium'],
                ['Bacteria',
                 'Firmicutes',
                 'Alicyclobacillaceae',
                 'Bacilli',
                 'Lactobacillales',
                 'Lactobacillales',
                 'Streptococcaceae',
                 'Streptococcus'],
                ['Bacteria',
                 'Actinobacteria',
                 'Actinobacteridae',
                 'Gordoniaceae',
                 'Corynebacteriaceae'],
                ['Bacteria',
                 'Firmicutes',
                 'Alicyclobacillaceae',
                 'Bacilli',
                 'Staphylococcaceae'],
                ['Bacteria', 'Cyanobacteria', 'Chloroplasts', 'vectors']])
        self.assertEqual(obs[0], exp[0])
        self.assertEqual(obs[1], exp[1])
        assert_almost_equal(obs[2], exp[2])
        self.assertEqual(obs[3], exp[3])

    def test_make_envs_dict(self):
        """ make_envs_dict should have the same abundance for each taxon
        as the matrix that made the dict"""
        envs = make_envs_dict(self.l19_data, self.l19_sample_names,
                              self.l19_taxon_names)
        for key in envs.keys():
            col_idx = self.l19_taxon_names.index(key)
            self.assertEqual(sum(envs[key].values()),
                             self.l19_data[:, col_idx].sum())

    def test_fields_to_dict(self):
        """fields_to_dict should make first field key, rest val"""
        test_data = \
            """0	R27DLI_4812	R27DLI_600	R27DLI_727	U1PLI_403	U1PLI_8969	U1PLI_9080	U1PLI_9526	W3Cecum_6642	W3Cecum_8992
1	U1PLI_7889
2	W3Cecum_4858
3	R27DLI_3243	R27DLI_4562	R27DLI_6828	R27DLI_9097	U1PLI_2780	U1PLI_67	U9PSI_10475	U9PSI_4341	W3Cecum_5191""".splitlines()  # output from cd-hit
        obs = fields_to_dict(test_data)
        exp = {
            '0': ['R27DLI_4812', 'R27DLI_600', 'R27DLI_727', 'U1PLI_403',
                  'U1PLI_8969', 'U1PLI_9080', 'U1PLI_9526', 'W3Cecum_6642', 'W3Cecum_8992'],
            '1': ['U1PLI_7889'],
            '2': ['W3Cecum_4858'],
            '3': ['R27DLI_3243', 'R27DLI_4562', 'R27DLI_6828', 'R27DLI_9097', 'U1PLI_2780', 'U1PLI_67', 'U9PSI_10475', 'U9PSI_4341', 'W3Cecum_5191']}
        self.assertEqual(obs, exp)

    def test_parse_qiime_parameters(self):
        """parse_qiime_parameters: functions with valid input """
        lines = ["#Don't edit this file!",
                 "pick_otus:similarity 0.94#this is not a comment...",
                 "pick_otus:otu_picking_method\tcdhit  # useful comment  ",
                 "align_seqs:verbose",
                 "assign_taxonomy:use_rdp\ttRuE # another great ## comment!",
                 "assign_taxonomy:something\tNone",
                 "",
                 "#some_script:fake_parameter\t99.0",
                 'summarize_taxa:md_identifier "Consensus Lineage"']
        actual = parse_qiime_parameters(lines)
        expected = {'pick_otus':
                    {'similarity': '0.94#this is not a comment...',
                     'otu_picking_method': 'cdhit'},
                    'assign_taxonomy':
                    {'use_rdp': None},
                    'summarize_taxa':
                    {'md_identifier': '"Consensus Lineage"'}}
        self.assertEqual(actual, expected)

        # default dict functions as expected -- looking up non-existant key
        # returns empty dict
        self.assertEqual(actual['some_other_script'], {})

    def test_parse_taxonomy(self):
        """ should parse taxonomy example, keeping otu id only"""
        example_tax = \
            """412 PC.635_647	Root;Bacteria;Firmicutes; "Clostridia";Clostridiales	0.930
319 PC.355_281	Root;Bacteria;Bacteroidetes   	0.970
353 PC.634_154	Root; Bacteria ; Bacteroidetes	0.830
17 PC.607_302	Root;Bacteria;Bacteroidetes	0.960
13 PC.481_1214	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales	0.870
338 PC.593_1314	Root; Bacteria 	0.990	42556	Additional fields ignored"""
        res = parse_taxonomy(example_tax.split('\n'))
        self.assertEqual(res['412'],
                         ["Root", "Bacteria", "Firmicutes", "\"Clostridia\"", "Clostridiales"])
        self.assertEqual(res['338'],
                         ["Root", "Bacteria"])

    def test_parse_taxonomy_to_otu_metadata(self):
        """parsing of taxonomy file to otu metadata format functions as expected
        """
        example_tax = \
            """412 PC.635_647	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales	0.930
319 PC.355_281	Root;Bacteria;Bacteroidetes	0.970
353 PC.634_154	Root;Bacteria;Bacteroidetes	0.830
17 PC.607_302	Root;Bacteria;Bacteroidetes	0.960
13	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales	0.870
338 PC.593_1314	Root;Bacteria	0.990"""
        actual = parse_taxonomy_to_otu_metadata(example_tax.split('\n'))
        expected = {
            '412': {'taxonomy': ['Root', 'Bacteria', 'Firmicutes',
                                 '"Clostridia"', 'Clostridiales'], 'score': 0.930},
            '319':
            {'taxonomy': ['Root',
                          'Bacteria',
                          'Bacteroidetes'],
             'score': 0.970},
            '353':
            {'taxonomy': ['Root',
                          'Bacteria',
                          'Bacteroidetes'],
             'score': 0.830},
            '17':
            {'taxonomy': ['Root',
                          'Bacteria',
                          'Bacteroidetes'],
             'score': 0.960},
            '13': {'taxonomy': ['Root', 'Bacteria', 'Firmicutes',
                                '"Clostridia"', 'Clostridiales'], 'score': 0.870},
            '338': {'taxonomy': ['Root', 'Bacteria'], 'score': 0.990}}
        self.assertEqual(actual, expected)

    def test_parse_taxonomy_to_otu_metadata_alt_labels(self):
        """parsing of taxonomy file to otu metadata format functions as expected
        """
        def f(v):
            return 1. + float(v)
        example_tax = \
            """412 PC.635_647	0.0
319 PC.355_281	0.970
353 PC.634_154	0.830
17 PC.607_302	0.960
13	0.870
338 PC.593_1314	0.990"""
        actual = parse_taxonomy_to_otu_metadata(
            example_tax.split('\n'),
            labels=['something'],
            process_fs=[f])
        expected = {'412': {'something': 1.0},
                    '319': {'something': 1.970},
                    '353': {'something': 1.830},
                    '17': {'something': 1.960},
                    '13': {'something': 1.870},
                    '338': {'something': 1.990}}
        self.assertEqual(actual, expected)

    def test_parse_taxonomy_to_otu_metadata_extra_fields_ignored(self):
        """parsing of taxonomy file to otu metadata format functions as expected
        """
        example_tax = \
            """412 PC.635_647	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales	0.930
319 PC.355_281	Root;Bacteria;Bacteroidetes	some text
353 PC.634_154	Root;Bacteria;Bacteroidetes	0.830
17 PC.607_302	Root;Bacteria;Bacteroidetes	0.960
13	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales	0.870
338 PC.593_1314	Root;Bacteria	0.990	42556	Additional fields ignored"""
        expected = {
            '412':
            {'taxonomy': ['Root',
                          'Bacteria',
                          'Firmicutes',
                          '"Clostridia"',
                          'Clostridiales']},
            '319': {'taxonomy': ['Root', 'Bacteria', 'Bacteroidetes']},
            '353': {'taxonomy': ['Root', 'Bacteria', 'Bacteroidetes']},
            '17': {'taxonomy': ['Root', 'Bacteria', 'Bacteroidetes']},
            '13':
            {'taxonomy': ['Root',
                          'Bacteria',
                          'Firmicutes',
                          '"Clostridia"',
                          'Clostridiales']},
            '338': {'taxonomy': ['Root', 'Bacteria']}}
        actual = parse_taxonomy_to_otu_metadata(
            example_tax.split('\n'),
            labels=['taxonomy'])
        self.assertEqual(actual, expected)

    def test_parse_taxonomy_to_otu_metadata_invalid_input(self):
        """parsing of taxonomy file to otu metadata format fails when too few functions
        """
        example_tax = \
            """412 PC.635_647	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales	0.930
319 PC.355_281	Root;Bacteria;Bacteroidetes	0.970
353 PC.634_154	Root;Bacteria;Bacteroidetes	0.830
17 PC.607_302	Root;Bacteria;Bacteroidetes	0.960
13	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales	0.870
338 PC.593_1314	Root;Bacteria	0.990"""
        self.assertRaises(
            ValueError,
            parse_taxonomy_to_otu_metadata,
            example_tax.split('\n'),
            labels=['taxonomy',
                    'score'],
            process_fs=[str])

    def test_parse_qiime_config_files(self):
        """ parse_qiime_config_files functions as expected """
        fake_file1 = ['key1\tval1', 'key2 val2']
        fake_file2 = ['key2\tval3']
        actual = parse_qiime_config_files([fake_file1, fake_file2])
        expected = {'key1': 'val1', 'key2': 'val3'}
        self.assertEqual(actual, expected)

        # looking up a non-existant value returns None
        self.assertEqual(actual['fake_key'], None)

        # empty dict on empty input
        self.assertEqual(parse_qiime_config_files([]), {})

        # test with an env variable - if it gets expanded
        # there won't be a $ in the output
        fake_file3 = ['key2\t$HOME', 'key3\thello $HOME']
        actual = parse_qiime_config_files([fake_file3])
        self.assertTrue('$' not in actual['key2'])
        self.assertTrue('$' not in actual['key3'])

    def test_parse_metadata_state_descriptions(self):
        """parse_metadata_state_descriptions should return correct states from string."""
        s = ''
        self.assertEqual(parse_metadata_state_descriptions(s), {})
        s = 'Study:Twin,Hand,Dog;BodySite:Palm,Stool'
        self.assertEqual(
            parse_metadata_state_descriptions(
                s), {'Study': set(['Twin', 'Hand', 'Dog']),
                     'BodySite': set(['Palm', 'Stool'])})

        # category names with colons i. e. ontology-derived
        s = 'Study:Twin,Hand,Dog;site:UBERON:feces,UBERON:ear canal;' +\
            'env_feature:ENVO:farm soil,ENVO:national park'
        self.assertEqual(parse_metadata_state_descriptions(s), {'Study':
                                                                set([
                                                                    'Twin', 'Hand', 'Dog']), 'site': set(['UBERON:feces',
                                                                                                          'UBERON:ear canal']), 'env_feature': set(['ENVO:farm soil',
                                                                                                                                                    'ENVO:national park'])})

        s = "Treatment:A,B,C;env_matter:ENVO:nitsol,ENVO:farm soil;env_biom:" +\
            "ENVO:Tropical dry (including Monsoon forests) and woodlands," +\
            "ENVO:Forest: including woodlands;country:GAZ:Persnickety Islands" +\
            ",St. Kitt's and Nevis"
        self.assertEqual(parse_metadata_state_descriptions(s), {"country":
                                                                set([
                                                                    "GAZ:Persnickety Islands", "St. Kitt's and Nevis"]),
                                                                "env_biom": set(["ENVO:Tropical dry (including Monsoon forests) " +
                                                                                 "and woodlands", "ENVO:Forest: including woodlands"]), "env_matter":
                                                                set([
                                                                    "ENVO:nitsol", "ENVO:farm soil"]), 'Treatment': set(["A", "B",
                                                                                                                         "C"])})

    def test_parse_illumina_line_barcode_in_header(self):
        """parse_illumina_line: handles barcode in header correctly """
        illumina_line0 = illumina_read1[0]
        illumina_line1 = illumina_read1[1]
        actual = parse_illumina_line(
            illumina_line0, barcode_length=6, rev_comp_barcode=True)
        expected = {
            'Full description': 'HWI-6X_9267:1:1:4:1699#ACCACCC/1',
            'Machine Name': 'HWI-6X_9267',
            'Channel Number': 1,
            'Tile Number': 1,
            'X Position': 4,
            'Y Position': 1699,
            'Barcode': 'GGTGGT',
            'Full Y Position Field': '1699#ACCACCC/1',
            'Sequence':
            'TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCC' +
            'CCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA',
            'Quality Score':
            'abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa' +
            'aaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB'}
        self.assertEqual(actual, expected)

        actual = parse_illumina_line(
            illumina_line0, barcode_length=6, rev_comp_barcode=False)
        expected['Barcode'] = 'ACCACC'

        actual = parse_illumina_line(
            illumina_line1, barcode_length=6, rev_comp_barcode=True)
        expected = {
            'Full description': 'HWI-6X_9267:1:1:4:390#ACCTCCC/1',
            'Machine Name': 'HWI-6X_9267',
            'Channel Number': 1,
            'Tile Number': 1,
            'X Position': 4,
            'Y Position': 390,
            'Barcode': 'GGAGGT',
            'Full Y Position Field': '390#ACCTCCC/1',
            'Sequence':
            'GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGG' +
            'GGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAA',
            'Quality Score':
            'aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa' +
            'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaBaaaaa'}
        self.assertEqual(actual, expected)

        actual = parse_illumina_line(
            illumina_line1, barcode_length=6, rev_comp_barcode=False)
        expected['Barcode'] = 'ACCTCC'

    def test_parse_illumina_line_barcode_in_sequence(self):
        """parse_illumina_line: handles barcode in sequence correctly """
        illumina_line0 = illumina_read3[0]
        actual = parse_illumina_line(
            illumina_line0, barcode_length=12,
            rev_comp_barcode=False, barcode_in_sequence=True)
        expected = {
            'Full description': 'HWI-EAS440_0386:1:23:19516:1031#0/1',
            'Machine Name': 'HWI-EAS440_0386',
            'Channel Number': 1,
            'Tile Number': 23,
            'X Position': 19516,
            'Y Position': 1031,
            'Barcode': 'ACAGCTAGCTTG',
            'Full Y Position Field': '1031#0/1',
            'Sequence':
            'TACGNAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATT'
            'GTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGGTGTCTTGAGT'
            'ACAGTAGAGGCAGGCGGAATTCGTGGGG',
            'Quality Score':
            'fffcGddd\_``_gggggggggggfgggggegggggcgggggggggggeeffafdcdfgbdggdbe]fbf'
            'dddddbdadadcddaf`abb`cVNRNUScaa``aOY]]]_[_BBBBBBBBBBB'
            'BBBBBBBBBBBBBBBBBBBBBBBBBBBBB'}
        self.assertEqual(actual, expected)

    def test_parse_qual_score(self):
        """qual_score should return dict of {id: qual_scores}"""
        scores = StringIO('>x\n5 10 5\n12\n>y\n30 40')
        self.assertItemsEqual(parse_qual_score(scores),
                         {'x': [5, 10, 5, 12], 'y': [30, 40]})

        # Check that a bad file, e.g. a fast raises Error
        bad_scores = StringIO('>x\nabcbd\n12\n>y\GATC')
        self.assertRaises(QiimeParseError, parse_qual_score, bad_scores)

    def test_parse_qual_scores(self):
        """qual_scores should return dict of {id:qual_scores}"""
        scores = StringIO('>x\n5 10 5\n12\n>y\n30 40')
        scores2 = StringIO('>a\n5 10 5\n12\n>b\n30 40')
        self.assertItemsEqual(parse_qual_scores([scores, scores2]),
                         {'x': [5, 10, 5, 12], 'y': [30, 40], 'a': [5, 10, 5, 12], 'b': [30, 40]})

    def test_MinimalQualParser(self):
        """MinimalQualParser should yield (id_, quals)"""
        scores = ['>x', '5 10 5', '12',
                  '>y', '30 40',
                  '>a', '5 10 5', '12',
                  '>b', '30 40']
        gen = list(MinimalQualParser(scores))
        self.assertItemsEqual(gen[0][1], [5, 10, 5, 12])
        self.assertItemsEqual(gen[1][1], [30, 40])
        self.assertItemsEqual(gen[2][1], [5, 10, 5, 12])
        self.assertItemsEqual(gen[3][1], [30, 40])

    def test_parse_trflp(self):
        """ should return a header and otu_table lists"""

        data = \
            """	Bin (10bp)	Bin (20bp)	Bin (30bp)	Bin (40 bp)
Samp-le 1	1000	2000	3000	4000
Sample 2		2000	3000	4000
Sample 3			3000	4000
Sample 4				4000
Sample 5	25			"""
        samples, otus, data = parse_trflp(data.split('\n'))

        samples_exp = [
            'Samp.le.1',
            'Sample.2',
            'Sample.3',
            'Sample.4',
            'Sample.5']
        otus_exp = ['Bin__10bp_', 'Bin__20bp_', 'Bin__30bp_', 'Bin__40_bp_']
        data_exp = array([[1000, 0, 0, 0, 25],
                          [2000, 2000, 0, 0, 0],
                          [3000, 3000, 3000, 0, 0],
                          [4000, 4000, 4000, 4000, 0]])

        self.assertEqual(samples, samples_exp)
        self.assertEqual(otus, otus_exp)
        assert_almost_equal(data, data_exp)

    def test_parse_trflp_headerless(self):
        """ should return a header and otu_table lists"""

        data = \
            """Samp-le 1	1000	2000	3000	4000
Sample 2		2000	3000	4000
Sample 3			3000	4000
Sample 4				4000
Sample_5__	25			"""
        samples, otus, data = parse_trflp(data.split('\n'))

        samples_exp = [
            'Samp.le.1',
            'Sample.2',
            'Sample.3',
            'Sample.4',
            'Sample.5..']
        otus_exp = ['Bin__0', 'Bin__1', 'Bin__2', 'Bin__3']
        data_exp = array([[1000, 0, 0, 0, 25],
                          [2000, 2000, 0, 0, 0],
                          [3000, 3000, 3000, 0, 0],
                          [4000, 4000, 4000, 4000, 0]])

        self.assertEqual(samples, samples_exp)
        self.assertEqual(otus, otus_exp)
        assert_almost_equal(data, data_exp)

    def test_parse_trflp_headerless_diff_row_len(self):
        """ should return a header and otu_table lists"""

        data = \
            """Samp-le 1	1000	2000
Sample 2		2000
Sample 3			3000
Sample 4				4000
Sample 5	25

"""
        samples, otus, data = parse_trflp(data.split('\n'))

        samples_exp = [
            'Samp.le.1',
            'Sample.2',
            'Sample.3',
            'Sample.4',
            'Sample.5']
        otus_exp = ['Bin__0', 'Bin__1', 'Bin__2', 'Bin__3']
        data_exp = array([[1000, 0, 0, 0, 25],
                          [2000, 2000, 0, 0, 0],
                          [0, 0, 3000, 0, 0],
                          [0, 0, 0, 4000, 0]])

        self.assertItemsEqual(samples, samples_exp)
        self.assertItemsEqual(otus, otus_exp)
        assert_almost_equal(data, data_exp)

    def test_parse_denoiser_mapping(self):
        """ parse_denoiser_mapping creates {} from denoiser mapping file
        """
        actual = parse_denoiser_mapping(self.denoiser_mapping1)
        expected = {'Read1': ['Read1', 'Read4', 'Read5 some comment'],
                    'Read2': ['Read2'],
                    'Read3': ['Read3', 'Read6']}
        self.assertDictEqual(actual, expected)

    def test_parse_otu_map(self):
        """ parse_otu_map functions as expected
        """
        otu_map_f = """otu1	s1_0	s2_1	s1_99
2	s1_9	s5_2 comment	s3_99	1_3	s1_75
otu3	s8_7	s2_5""".split('\n')
        expected_map = {(0, 0): 2, (0, 1): 1,
                        (1, 0): 2, (1, 2): 1, (1, 3): 1, (1, 4): 1,
                        (2, 5): 1, (2, 1): 1}
        expected_sids = ['s1', 's2', 's5', 's3', '1', 's8']
        expected_oids = ['otu1', '2', 'otu3']
        actual = parse_otu_map(otu_map_f)
        self.assertDictEqual(actual[0], expected_map)
        self.assertItemsEqual(actual[1], expected_sids)
        self.assertItemsEqual(actual[2], expected_oids)

    def test_parse_otu_map_w_excludes(self):
        """ parse_otu_map functions as expected when excluding otu ids
        """
        otu_map_f = """otu1	s1_0	s2_1	s1_99
2	s1_9	s5_2 comment	s3_99	1_3	s1_75
otu3	s8_7	s2_5""".split('\n')
        excludes = ['otu1', '2']
        expected_map = {(0, 0): 1, (0, 1): 1}
        expected_sids = ['s8', 's2']
        expected_oids = ['otu3']
        actual = parse_otu_map(otu_map_f, excludes)
        self.assertEqual(actual[0], expected_map)
        self.assertEqual(actual[1], expected_sids)
        self.assertEqual(actual[2], expected_oids)

    def test_parse_sample_id_map(self):
        """Test parsing a sample id map functions correctly."""
        sample_id_map = ['\t\t\n', '', ' ', '\n', 'S1\ta',
                         'S2\tb', '\n \t', 'T1\ta', 'T2\tb']
        exp = {'S1': 'a', 'S2': 'b', 'T1': 'a', 'T2': 'b'}
        obs = parse_sample_id_map(sample_id_map)
        self.assertEqual(obs, exp)

    def test_parse_sample_id_map_repeat_sample_ids(self):
        """Test parsing a sample id map with non-unique first column fails."""
        sample_id_map = ['\t\t\n', '', ' ', '\n', 'S1\ta',
                         'S2\tb', '\n \t', 'S1\tc']
        self.assertRaises(ValueError, parse_sample_id_map,
                          sample_id_map)

    def test_parse_sample_id_map_many_to_one_mapping(self):
        """Test parsing a sample id map with many-to-one mapping fails."""
        sample_id_map = ['S1\ta', 'T1\ta', 'S2\ta']
        self.assertRaises(ValueError, parse_sample_id_map,
                          sample_id_map)


illumina_read1 = """HWI-6X_9267:1:1:4:1699#ACCACCC/1:TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA:abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB
HWI-6X_9267:1:1:4:390#ACCTCCC/1:GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAA:aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaBaaaaa""".split('\n')

illumina_read2 = """HWI-6X_9267:1:1:4:1699#ACCACCC/2:TTTTAAAAAAAAGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTAAAAAAAAACCCCCCCGGGGGGGGTTTTTTTAATTATTC:aaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbcccccccccccccccccBcccccccccccccccc```````BBBB
HWI-6X_9267:1:1:4:390#ACCTCCC/2:ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG:aaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbb""".split('\n')

illumina_read3 = """HWI-EAS440_0386:1:23:19516:1031#0/1:ACAGCTAGCTTGTACGNAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGGTGTCTTGAGTACAGTAGAGGCAGGCGGAATTCGTGGGG:gggggggeggcffffcGddd\_``_gggggggggggfgggggegggggcgggggggggggeeffafdcdfgbdggdbe]fbfdddddbdadadcddaf`abb`cVNRNUScaa``aOY]]]_[_BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
HWI-EAS440_0386:1:23:19660:1034#0/1:CATATCGCAGTTTACGNAAGGTCCGAGCGTTGTCCGGAATCATTGGGCGTAAAGGGTACGTAGGCGGGTAAGCAAGTTAGAAGTGAAATCCTATAGCTCAACTATAGTAAGCTTTTAAAACTGCTCATCTTGAGGTATGGAAGGGAAAGTGGAATTCCTAGTTA:fhhghhhhfhghhhhcHdcddccddhhhhhhfhhhhghhhdghhhhhhhhhhfhhhdhghgghhhhhhbfdfdbagdgdgfffafa]dad_acdabZcaabad[a__^_`cbddefb_cd^]_L\]U_]^aaZ___]bBBBBBBBBBBBBBBBBBBBBBBBBBB""".split('\n')

legacy_otu_table1 = """# some comment goes here
#OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium

1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

otu_table1 = """# Some comment




OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
# some other comment
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
# comments
#    everywhere!
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

otu_table1_floats = """# Some comment




OTU ID	Fing	Key	NA	Consensus Lineage
0	19111.0	44536.0	42.0	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
# some other comment
1	1216.0	3500.0	6.0	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803.0	1184.0	2.0	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
# comments
#    everywhere!
3	1722.1	4903.2	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589.6	2074.4	34.5	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""


otu_table_without_leading_comment = '#OTU ID\tlet-7i\tmiR-7\tmiR-17n\tmiR-18a\tmiR-19a\tmiR-22\tmiR-25\tmiR-26a\nA2M\t-0.2\t0.03680505\t0.205\t0.23\t0.66\t0.08\t-0.373\t0.26\nAAAS\t-0.09\t-0.25\t0.274\t0.15\t0.12\t0.29\t0.029\t-0.114845199\nAACS\t0.33\t0.19\t0.27\t0.28\t0.19\t0.25\t0.089\t0.14\nAADACL1\t0.49\t-0.92\t-0.723\t-0.23\t0.08\t0.49\t-0.386\t-0.64'


taxa_summary1 = """#Full OTU Counts
Taxon	Even1	Even2	Even3
Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae	0.0880247251673	0.0721968465746	0.081371761759
Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia	0.192137761955	0.191095101593	0.188504131885
Bacteria;Firmicutes;Bacilli;Lactobacillales	0.0264895739603	0.0259942669171	0.0318460745596
# some comment
Bacteria;Firmicutes;Clostridia;Clostridiales	0.491800007824	0.526186212556	0.49911159984
Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales	0.0311411916592	0.0184083913576	0.0282325481054
Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales	0.166137214246	0.163087129528	0.168923372865
No blast hit;Other	0.00426952518811	0.00303205147361	0.0020105109874"""

taxa_summary1_expected = (['Even1', 'Even2', 'Even3'],
                          ['Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae',
                           'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
                           'Bacteria;Firmicutes;Bacilli;Lactobacillales',
                           'Bacteria;Firmicutes;Clostridia;Clostridiales',
                           'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
                           'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
                           'No blast hit;Other'],
                          array(
                              [[0.0880247251673, 0.0721968465746, 0.081371761759],
                               [0.192137761955,
                                0.191095101593,
                                0.188504131885],
                               [0.0264895739603,
                                0.0259942669171,
                                0.0318460745596],
                               [0.491800007824, 0.526186212556, 0.49911159984],
                               [0.0311411916592,
                                0.0184083913576,
                                0.0282325481054],
                               [0.166137214246,
                                0.163087129528,
                                0.168923372865],
                               [0.00426952518811, 0.00303205147361, 0.0020105109874]]))

expected_lineages1 = [
    ['Bacteria',
     'Actinobacteria',
     'Actinobacteridae',
     'Propionibacterineae',
     'Propionibacterium'],
    ['Bacteria',
     'Firmicutes',
     'Alicyclobacillaceae',
     'Bacilli',
     'Lactobacillales',
     'Lactobacillales',
     'Streptococcaceae',
     'Streptococcus'],
    ['Bacteria',
     'Actinobacteria',
     'Actinobacteridae',
     'Gordoniaceae',
     'Corynebacteriaceae'],
    ['Bacteria',
     'Firmicutes',
     'Alicyclobacillaceae',
     'Bacilli',
     'Staphylococcaceae'],
    ['Bacteria', 'Cyanobacteria', 'Chloroplasts', 'vectors']]

denoiser_mapping1 = """Read1:\tRead4\tRead5 some comment
Read2:
Read3:\tRead6"""

sam_data1 = """@SQ	SN:s1_1	LN:66
@SQ	SN:s2_2	LN:1131
@SQ	SN:s1_3	LN:348
@SQ	SN:s1_4	LN:348
@SQ	SN:s1_5	LN:1131
@SQ	SN:s1_6	LN:1132
s1_1	0	s1_1	2	136	1S65M	*	0	0	atgaaacgcattagcaccaccattaccaccaccatcaccattaccacaggtaacggtgcgggctga	*	AS:i:65	XS:i:0	XF:i:3	XE:i:2	NM:i:0
s2_2	0	s1_5	1	0	1131M	*	0	0	atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaaatcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgacaaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcgcaaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggcggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgatatttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataacatggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgactctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagacttgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacagcagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgtcatggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtggacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccggcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggcaacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatcgaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggtaagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgatttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctgctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctcaaagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa	*	AS:i:1131	XS:i:1131	XF:i:0	XE:i:24	NM:i:0
s1_3	0	s1_3	1	9	348M	*	0	0	atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagttacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtctttctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtctttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagccctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgcaaagccctggcccgcgattccttaagcctgcttgcctacgtcaaataa	*	AS:i:348	XS:i:336	XF:i:0	XE:i:9	NM:i:0
s1_4	0	s1_4	1	9	348M	*	0	0	atgaagaaaattttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagttacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtctttctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtctttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagccctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgcaaagccctggcccgcgattccttaagcctgcttgcctacgtcaaatcc	*	AS:i:348	XS:i:336	XF:i:0	XE:i:9	NM:i:0
s1_5	0	s2_2	1	0	1131M	*	0	0	atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaaatcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgacaaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcgcaaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggcggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgatatttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataacatggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgactctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagacttgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacagcagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgtcatggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtggacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccggcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggcaacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatcgaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggtaagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgatttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctgctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctcaaagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa	*	AS:i:1131	XS:i:1131	XF:i:0	XE:i:24	NM:i:0
s1_6	0	s1_6	1	1	1132M	*	0	0	aatgactaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaaatcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgacaaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcgcaaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggcggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgatatttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataacatggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgactctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagacttgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacagcagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgtcatggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtggacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccggcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggcaacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatcgaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggtaagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgatttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctgctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctcaaagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa	*	AS:i:1132	XS:i:1128	XF:i:0	XE:i:24	NM:i:0
"""

sam1_expected = [
    ["s1_1",
     "0",
     "s1_1",
     "2",
     "136",
     "1S65M",
     "*",
     "0",
     "0",
     "atgaaacgcattagcaccaccattaccaccaccatcaccattaccacaggtaacggtgcgggctga",
     "*",
     "AS:i:65",
     "XS:i:0",
     "XF:i:3",
     "XE:i:2",
     "NM:i:0"],
    ["s2_2",
     "0",
     "s1_5",
     "1",
     "0",
     "1131M",
     "*",
     "0",
     "0",
     "atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaaatcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgacaaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcgcaaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggcggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgatatttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataacatggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgactctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagacttgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacagcagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgtcatggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtggacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccggcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggcaacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatcgaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggtaagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgatttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctgctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctcaaagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa",
     "*",
     "AS:i:1131",
     "XS:i:1131",
     "XF:i:0",
     "XE:i:24",
     "NM:i:0"],
    ["s1_3",
     "0",
     "s1_3",
     "1",
     "9",
     "348M",
     "*",
     "0",
     "0",
     "atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagttacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtctttctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtctttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagccctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgcaaagccctggcccgcgattccttaagcctgcttgcctacgtcaaataa",
     "*",
     "AS:i:348",
     "XS:i:336",
     "XF:i:0",
     "XE:i:9",
     "NM:i:0"],
    ["s1_4",
     "0",
     "s1_4",
     "1",
     "9",
     "348M",
     "*",
     "0",
     "0",
     "atgaagaaaattttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagttacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtctttctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtctttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagccctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgcaaagccctggcccgcgattccttaagcctgcttgcctacgtcaaatcc",
     "*",
     "AS:i:348",
     "XS:i:336",
     "XF:i:0",
     "XE:i:9",
     "NM:i:0"],
    ["s1_5",
     "0",
     "s2_2",
     "1",
     "0",
     "1131M",
     "*",
     "0",
     "0",
     "atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaaatcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgacaaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcgcaaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggcggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgatatttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataacatggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgactctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagacttgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacagcagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgtcatggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtggacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccggcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggcaacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatcgaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggtaagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgatttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctgctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctcaaagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa",
     "*",
     "AS:i:1131",
     "XS:i:1131",
     "XF:i:0",
     "XE:i:24",
     "NM:i:0"],
    ["s1_6", "0", "s1_6", "1", "1", "1132M", "*", "0", "0", "aatgactaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaaatcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgacaaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcgcaaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggcggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgatatttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataacatggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgactctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagacttgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacagcagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgtcatggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtggacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccggcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggcaacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatcgaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggtaagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgatttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctgctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctcaaagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa", "*", "AS:i:1132", "XS:i:1128", "XF:i:0", "XE:i:24", "NM:i:0"]]

COORDS_NO_VECTORS = """A\t0.11\t0.09\t0.23
B\t0.03\t0.07\t-0.26
C\t0.12\t0.06\t-0.32
eigvals\t4.94\t1.79\t1.50
% variation explained\t14.3\t5.2\t4.3"""

COORDS_NO_EIGENVALS = """pc vector number\t1\t2\t3
A\t0.11\t0.09\t0.23
B\t0.03\t0.07\t-0.26
C\t0.12\t0.06\t-0.32
foo\t4.94\t1.79\t1.50
% variation explained\t14.3\t5.2\t4.3"""

COORDS_NO_PCNTS = """pc vector number\t1\t2\t3
A\t0.11\t0.09\t0.23
B\t0.03\t0.07\t-0.26
C\t0.12\t0.06\t-0.32
eigvals\t4.94\t1.79\t1.50"""

if __name__ == '__main__':
    main()
