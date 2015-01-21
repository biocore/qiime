#!/usr/bin/env python

from os import close
from tempfile import mkstemp

from unittest import TestCase, main
from biom import load_table

from skbio.util import remove_files
from qiime.util import load_qiime_config
from qiime.differential_abundance import DA_fitZIG, DA_DESeq2, check_mapping_file_category


__author__ = "Sophie Weiss"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Sophie Weiss"]
__license__ = "GPL"
__version__ = "1.9.0-rc1"
__maintainer__ = "Sophie Weiss"
__email__ = "sophie.sjw@gmail.com"

"""Contains tests for the output of DESeq2 negative binomial and metagenomeSeq fitZIG OTU differential abundance testing.  modeled after test_supervised_learning.py"""


def is_float(input_string):
    """True if string can be cast as a float"""
    try:
        float(input_string)
        return True
    except ValueError:
        return False

class RDifferentialAbundanceTests(TestCase):

    """Tests of the RDifferentialAbundanceTest class"""

    def setUp(self):
        self.qiime_config = load_qiime_config()
        self.tmp_dir = self.qiime_config['temp_dir'] or '/tmp/'

        # Temporary input file
        fd, self.tmp_otu_fp = mkstemp(dir=self.tmp_dir,
                                        prefix='R_test_otu_table_',
                                           suffix='.biom')
        close(fd)
        seq_file = open(self.tmp_otu_fp, 'w')
        seq_file.write(test_otu_table)
        seq_file.close()

        fd, self.tmp_map_fp = mkstemp(dir=self.tmp_dir,
                                        prefix='R_test_map_',
                                           suffix='.txt')
        close(fd)
        seq_file = open(self.tmp_map_fp, 'w')
        seq_file.write(test_map)
        seq_file.close()

        fd, self.tmp_otu_fp_fitZIG_out = mkstemp(dir=self.tmp_dir,
                                                prefix='R_test_otu_table_fitZIG_out_',
                                                suffix='.txt')
        fd, self.tmp_otu_fp_DESeq2_out = mkstemp(dir=self.tmp_dir,
                                                prefix='R_test_otu_table_DESeq2_out_',
                                                suffix='.txt')

        self.files_to_remove = \
            [self.tmp_otu_fp, self.tmp_map_fp, self.tmp_otu_fp_fitZIG_out, self.tmp_otu_fp_DESeq2_out]

        DA_fitZIG(self.tmp_otu_fp, self.tmp_otu_fp_fitZIG_out, self.tmp_map_fp, 'Individual', 'S1', 'S2')
        DA_DESeq2(self.tmp_otu_fp, self.tmp_otu_fp_DESeq2_out, self.tmp_map_fp, 'Individual', 'S1', 'S2', DESeq2_diagnostic_plots=False)


    def test_check_mapping_file_category(self):
        z = load_table(self.tmp_otu_fp)

        with self.assertRaises(ValueError):
            check_mapping_file_category(z, self.tmp_map_fp, 'Q', 'S1', 'S2')

        with self.assertRaises(ValueError):
            check_mapping_file_category(z, self.tmp_map_fp, 'Individual', 'dog', 'S2')

        with self.assertRaises(ValueError):
            check_mapping_file_category(z, self.tmp_map_fp, 'Individual', 'S1', 'dog')

        with self.assertRaises(ValueError):
            check_mapping_file_category(z, self.tmp_map_fp, 'Individual', 'S1', 'S1')


    def test_metagenomeSeq_fitZIG_format(self):
        zig = open(self.tmp_otu_fp_fitZIG_out).readlines()
        #test header format
        exp = 'OTU\t+samples in group 0\t+samples in group 1\tcounts in group 0\tcounts in group 1\toddsRatio\tlower\tupper\tfisherP\tfisherAdjP\t(Intercept)\tMGS_categoryS2\tMGS_categoryS3\tscalingFactor\tpvalues\tadjPvalues\ttaxonomy\n'
        self.assertEqual(zig[0], exp)

        #test that fitZIG returns 25 features, it only returns the most important
        num_features_returned = len(zig) - 1
        self.assertEqual(num_features_returned, 25)

        # ensure that each line has two elements, and that the first one
        # is the name of one of the OTUs, the second is a float
        for line in zig[1:]:
            words = line.strip().split('\t')
            line_length = len(words)
            self.assertEqual(line_length, 17)
            self.assertEqual(words[0] in test_OTU_IDs, True)
            self.assertEqual(is_float(words[1]), True)

        #test first five significant OTUs
        exp = ['1487\t3\t13\t3\t54\t0\t0\t0.0814754160126085\t1.03417385089143e-07\t2.58543462722858e-06\t0.732914457135622\t2.35091578797182\t-0.358285072078857\t-28.0403903601361\t2.59620832179624e-05\t0.000649052080449059\tAy; other\n',
               '193\t2\t10\t4\t25\t0.0308298124876527\t0.00223706309126031\t0.229093863501484\t3.241507131278e-05\t0.000135062797136583\t0.980060835463344\t2.05824256013173\t-0.393172554812399\t-43.9276179169663\t6.26549562011497e-05\t0.000783186952514372\tAe; other\n',
               '979\t0\t10\t0\t20\t0\t0\t0.0956087104342576\t6.04991702771485e-07\t5.04159752309571e-06\t0.0855575301280384\t1.4678274119766\t0.00901541795158288\t-7.07862887279825\t0.000648824231851767\t0.00540686859876472\tAp; other\n',
               '1314\t4\t0\t7\t0\tInf\t0.34924924266116\tInf\t0.277924541082436\t0.347405676353045\t0.412980602558517\t-1.25800421238268\t-0.711956542078996\t22.3778549863952\t0.00315666606837298\t0.0197291629273311\tAv; other\n',
               '1351\t0\t6\t0\t12\t0\t0\t0.320469784786932\t0.000621585760904647\t0.00194245550282702\t0.106944586260241\t1.21145998904839\t0.0112690273007214\t-8.84809361558598\t0.00579444501011918\t0.0289722250505959\tAw; other\n']
        for a, e in zip(zig[1:6],exp):
            af = map(str,a.split('\t'))
            ef = map(str,e.split('\t'))
            self.assertEqual(len(af), len(ef))
            at = af.pop()
            et = ef.pop()
            self.assertEqual(at, et)
            af = map(float,af)
            ef = map(float,ef)
            for af_e, ef_e in zip(af, ef):
                self.assertAlmostEqual(af_e, ef_e)

    def test_DESeq2_nbinom_format(self):
        nbinom = open(self.tmp_otu_fp_DESeq2_out).readlines()
        #test header format
        exp = 'OTU\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\ttaxonomy\n'
        self.assertEqual(nbinom[0], exp)

        #test that all features in input biom table returned
        num_features_returned = len(nbinom) - 1
        self.assertEqual(num_features_returned, 50)

        # ensure that each line has seven elements, and that the first one
        # is the name of one of the OTUs, the second is a float
        for line in nbinom[1:]:
            words = line.strip().split('\t')
            line_length = len(words)
            self.assertEqual(line_length, 8)
            self.assertEqual(words[0] in test_OTU_IDs, True)
            self.assertEqual(is_float(words[1]), True)

        #test first five significant OTUs
        exp=['1848\t119.050214384222\t2.9560510780622\t0.27609765938806\t10.706541607828\t9.48387951847218e-27\t2.84516385554165e-25\tBc; other\n',
             '2096\t4.57782752174272\t-2.51043602053873\t0.534006762361075\t-4.7011315164606\t2.58723804969022e-06\t3.88085707453532e-05\tBg; other\n',
             '2366\t3.63266497414694\t-2.21660628894061\t0.514240606913112\t-4.31044584799803\t1.62925724104894e-05\t0.000162925724104894\tBm; other\n',
             '1088\t2.80571100244187\t-1.9610730170819\t0.485853003903793\t-4.03635050380428\t5.42890917069973e-05\t0.00040716818780248\tAs; other\n',
             '3006\t2.62956013500845\t-1.73771036728344\t0.489095213973481\t-3.55290814065839\t0.000380997482497004\t0.00228598489498202\tBu; other\n']
        for a, e in zip(nbinom[1:6],exp):
            af = map(str,a.split('\t'))
            ef = map(str,e.split('\t'))
            self.assertEqual(len(af), len(ef))
            at = af.pop()
            et = ef.pop()
            self.assertEqual(at, et)
            af = map(float,af)
            ef = map(float,ef)
            for af_e, ef_e in zip(af, ef):
                self.assertAlmostEqual(af_e, ef_e)

    def tearDown(self):
        remove_files(set(self.files_to_remove))


test_OTU_IDs = ['88',
                '131',
                '144',
                '158',
                '193',
                '225',
                '260',
                '588',
                '634',
                '721',
                '821',
                '843',
                '883',
                '891',
                '976',
                '979',
                '983',
                '1035',
                '1088',
                '1156',
                '1287',
                '1314',
                '1351',
                '1373',
                '1487',
                '1582',
                '1591',
                '1784',
                '1848',
                '1886',
                '2007',
                '2059',
                '2096',
                '2187',
                '2218',
                '2270',
                '2328',
                '2360',
                '2366',
                '2407',
                '2519',
                '2526',
                '2810',
                '2915',
                '2932',
                '2956',
                '3006',
                '3060',
                '3108',
                '3127']

test_otu_table = """{"rows": [{"id": "88", "metadata": {"taxonomy": ["Aa", "other"]}}, {"id": "131", "metadata": {"taxonomy": ["Ab", "other"]}}, {"id": "144", "metadata": {"taxonomy": ["Ac", "other"]}}, {"id": "158", "metadata": {"taxonomy": ["Ad", "other"]}}, {"id": "193", "metadata": {"taxonomy": ["Ae", "other"]}}, {"id": "225", "metadata": {"taxonomy": ["Af", "other"]}}, {"id": "260", "metadata": {"taxonomy": ["Ag", "other"]}}, {"id": "588", "metadata": {"taxonomy": ["Ah", "other"]}}, {"id": "634", "metadata": {"taxonomy": ["Ai", "other"]}}, {"id": "721", "metadata": {"taxonomy": ["Aj", "other"]}}, {"id": "821", "metadata": {"taxonomy": ["Ak", "other"]}}, {"id": "843", "metadata": {"taxonomy": ["Al", "other"]}}, {"id": "883", "metadata": {"taxonomy": ["Am", "other"]}}, {"id": "891", "metadata": {"taxonomy": ["An", "other"]}}, {"id": "976", "metadata": {"taxonomy": ["Ao", "other"]}}, {"id": "979", "metadata": {"taxonomy": ["Ap", "other"]}}, {"id": "983", "metadata": {"taxonomy": ["Aq", "other"]}}, {"id": "1035", "metadata": {"taxonomy": ["Ar", "other"]}}, {"id": "1088", "metadata": {"taxonomy": ["As", "other"]}}, {"id": "1156", "metadata": {"taxonomy": ["At", "other"]}}, {"id": "1287", "metadata": {"taxonomy": ["Au", "other"]}}, {"id": "1314", "metadata": {"taxonomy": ["Av", "other"]}}, {"id": "1351", "metadata": {"taxonomy": ["Aw", "other"]}}, {"id": "1373", "metadata": {"taxonomy": ["Ax", "other"]}}, {"id": "1487", "metadata": {"taxonomy": ["Ay", "other"]}}, {"id": "1582", "metadata": {"taxonomy": ["Az", "other"]}}, {"id": "1591", "metadata": {"taxonomy": ["Ba", "other"]}}, {"id": "1784", "metadata": {"taxonomy": ["Bb", "other"]}}, {"id": "1848", "metadata": {"taxonomy": ["Bc", "other"]}}, {"id": "1886", "metadata": {"taxonomy": ["Bd", "other"]}}, {"id": "2007", "metadata": {"taxonomy": ["Be", "other"]}}, {"id": "2059", "metadata": {"taxonomy": ["Bf", "other"]}}, {"id": "2096", "metadata": {"taxonomy": ["Bg", "other"]}}, {"id": "2187", "metadata": {"taxonomy": ["Bh", "other"]}}, {"id": "2218", "metadata": {"taxonomy": ["Bi", "other"]}}, {"id": "2270", "metadata": {"taxonomy": ["Bj", "other"]}}, {"id": "2328", "metadata": {"taxonomy": ["Bk", "other"]}}, {"id": "2360", "metadata": {"taxonomy": ["Bl", "other"]}}, {"id": "2366", "metadata": {"taxonomy": ["Bm", "other"]}}, {"id": "2407", "metadata": {"taxonomy": ["Bn", "other"]}}, {"id": "2519", "metadata": {"taxonomy": ["Bo", "other"]}}, {"id": "2526", "metadata": {"taxonomy": ["Bp", "other"]}}, {"id": "2810", "metadata": {"taxonomy": ["Bq", "other"]}}, {"id": "2915", "metadata": {"taxonomy": ["Br", "other"]}}, {"id": "2932", "metadata": {"taxonomy": ["Bs", "other"]}}, {"id": "2956", "metadata": {"taxonomy": ["Bt", "other"]}}, {"id": "3006", "metadata": {"taxonomy": ["Bu", "other"]}}, {"id": "3060", "metadata": {"taxonomy": ["Bv", "other"]}}, {"id": "3108", "metadata": {"taxonomy": ["Bw", "other"]}}, {"id": "3127", "metadata": {"taxonomy": ["Bx", "other"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 13, 1.0], [0, 14, 4.0], [0, 15, 2.0], [0, 17, 3.0], [0, 18, 4.0], [0, 19, 3.0], [0, 22, 1.0], [0, 23, 5.0], [0, 24, 3.0], [0, 25, 3.0], [0, 26, 1.0], [0, 27, 1.0], [0, 28, 1.0], [0, 29, 1.0], [0, 30, 2.0], [0, 31, 1.0], [0, 35, 1.0], [0, 37, 1.0], [1, 9, 1.0], [1, 12, 1.0], [1, 13, 1.0], [1, 14, 1.0], [1, 17, 1.0], [1, 19, 1.0], [1, 20, 1.0], [1, 22, 1.0], [2, 12, 2.0], [2, 13, 1.0], [2, 14, 4.0], [2, 16, 2.0], [2, 17, 4.0], [2, 18, 6.0], [2, 19, 7.0], [2, 23, 2.0], [2, 24, 4.0], [2, 37, 1.0], [3, 12, 2.0], [3, 13, 1.0], [3, 14, 3.0], [3, 15, 5.0], [3, 17, 1.0], [3, 18, 7.0], [3, 19, 6.0], [3, 20, 1.0], [3, 21, 1.0], [3, 23, 1.0], [3, 24, 3.0], [3, 26, 1.0], [3, 27, 1.0], [3, 35, 1.0], [4, 0, 2.0], [4, 4, 2.0], [4, 12, 3.0], [4, 13, 8.0], [4, 15, 1.0], [4, 16, 2.0], [4, 17, 2.0], [4, 18, 1.0], [4, 20, 1.0], [4, 21, 3.0], [4, 22, 3.0], [4, 24, 1.0], [5, 4, 1.0], [5, 8, 3.0], [5, 11, 2.0], [5, 12, 1.0], [5, 13, 1.0], [5, 21, 1.0], [5, 23, 1.0], [5, 26, 2.0], [5, 27, 2.0], [5, 29, 3.0], [5, 31, 2.0], [5, 34, 1.0], [5, 36, 3.0], [5, 37, 5.0], [6, 0, 1.0], [6, 3, 1.0], [6, 5, 2.0], [6, 7, 1.0], [6, 9, 1.0], [6, 37, 1.0], [7, 2, 1.0], [7, 20, 1.0], [7, 23, 1.0], [7, 25, 1.0], [7, 26, 1.0], [7, 29, 1.0], [7, 31, 1.0], [7, 32, 1.0], [7, 37, 1.0], [8, 7, 1.0], [8, 12, 3.0], [8, 15, 2.0], [8, 17, 1.0], [8, 18, 1.0], [8, 20, 2.0], [8, 21, 2.0], [8, 23, 2.0], [8, 28, 1.0], [8, 29, 1.0], [8, 32, 1.0], [8, 34, 1.0], [8, 36, 1.0], [9, 13, 1.0], [9, 14, 2.0], [9, 16, 1.0], [9, 18, 1.0], [9, 22, 1.0], [10, 0, 5.0], [10, 3, 1.0], [10, 4, 1.0], [10, 5, 1.0], [10, 6, 1.0], [10, 7, 2.0], [10, 9, 3.0], [10, 10, 1.0], [10, 11, 1.0], [10, 12, 6.0], [10, 13, 10.0], [10, 14, 4.0], [10, 15, 5.0], [10, 16, 2.0], [10, 17, 1.0], [10, 18, 6.0], [10, 19, 2.0], [10, 20, 3.0], [10, 21, 3.0], [10, 22, 10.0], [10, 23, 8.0], [10, 24, 3.0], [10, 27, 1.0], [10, 31, 2.0], [10, 32, 1.0], [10, 33, 1.0], [10, 34, 2.0], [10, 36, 1.0], [11, 0, 1.0], [11, 2, 1.0], [11, 3, 3.0], [11, 10, 2.0], [11, 12, 2.0], [11, 13, 4.0], [11, 15, 2.0], [11, 16, 1.0], [11, 17, 3.0], [11, 18, 4.0], [11, 20, 1.0], [11, 21, 2.0], [11, 23, 2.0], [11, 24, 3.0], [11, 29, 2.0], [12, 0, 4.0], [12, 1, 7.0], [12, 2, 3.0], [12, 3, 1.0], [12, 4, 1.0], [12, 5, 3.0], [12, 6, 5.0], [12, 7, 1.0], [12, 8, 1.0], [12, 9, 2.0], [12, 11, 5.0], [12, 13, 2.0], [12, 14, 3.0], [12, 15, 1.0], [12, 17, 1.0], [12, 18, 1.0], [12, 19, 3.0], [12, 20, 14.0], [12, 21, 10.0], [12, 22, 6.0], [12, 23, 3.0], [12, 24, 2.0], [12, 30, 1.0], [12, 34, 3.0], [12, 37, 1.0], [13, 1, 1.0], [13, 2, 5.0], [13, 6, 3.0], [13, 8, 1.0], [13, 11, 1.0], [13, 12, 5.0], [13, 13, 3.0], [13, 14, 10.0], [13, 15, 3.0], [13, 16, 5.0], [13, 17, 12.0], [13, 18, 10.0], [13, 19, 8.0], [13, 20, 2.0], [13, 22, 3.0], [13, 23, 3.0], [13, 24, 9.0], [13, 26, 2.0], [13, 27, 3.0], [13, 28, 8.0], [13, 29, 2.0], [13, 30, 8.0], [13, 31, 2.0], [13, 32, 1.0], [13, 33, 1.0], [13, 34, 1.0], [13, 35, 3.0], [13, 36, 4.0], [13, 37, 2.0], [14, 1, 4.0], [14, 4, 8.0], [14, 9, 1.0], [14, 11, 1.0], [14, 12, 4.0], [14, 13, 16.0], [14, 14, 2.0], [14, 15, 4.0], [14, 16, 1.0], [14, 18, 1.0], [14, 20, 1.0], [14, 23, 3.0], [14, 25, 7.0], [14, 26, 2.0], [14, 27, 9.0], [14, 28, 2.0], [14, 29, 3.0], [14, 30, 5.0], [14, 31, 17.0], [14, 32, 6.0], [14, 33, 10.0], [14, 34, 4.0], [14, 35, 7.0], [14, 36, 12.0], [14, 37, 6.0], [15, 12, 2.0], [15, 13, 2.0], [15, 14, 2.0], [15, 15, 2.0], [15, 16, 6.0], [15, 17, 1.0], [15, 19, 1.0], [15, 20, 1.0], [15, 23, 2.0], [15, 24, 1.0], [16, 0, 2.0], [16, 1, 1.0], [16, 5, 2.0], [16, 6, 1.0], [16, 9, 1.0], [16, 10, 1.0], [16, 12, 7.0], [16, 14, 3.0], [16, 20, 1.0], [16, 21, 1.0], [16, 23, 3.0], [16, 29, 1.0], [17, 0, 7.0], [17, 1, 10.0], [17, 2, 4.0], [17, 3, 3.0], [17, 4, 3.0], [17, 5, 19.0], [17, 6, 9.0], [17, 7, 11.0], [17, 8, 10.0], [17, 9, 5.0], [17, 10, 4.0], [17, 11, 13.0], [17, 12, 7.0], [17, 13, 18.0], [17, 14, 7.0], [17, 15, 16.0], [17, 16, 6.0], [17, 17, 21.0], [17, 18, 26.0], [17, 19, 16.0], [17, 20, 8.0], [17, 21, 23.0], [17, 22, 17.0], [17, 23, 9.0], [17, 24, 30.0], [17, 25, 4.0], [17, 26, 7.0], [17, 27, 7.0], [17, 28, 6.0], [17, 29, 6.0], [17, 30, 3.0], [17, 31, 4.0], [17, 32, 2.0], [17, 34, 2.0], [17, 35, 10.0], [17, 36, 2.0], [17, 37, 3.0], [18, 12, 2.0], [18, 13, 1.0], [18, 14, 4.0], [18, 15, 1.0], [18, 16, 9.0], [18, 17, 12.0], [18, 18, 7.0], [18, 19, 9.0], [18, 20, 3.0], [18, 21, 2.0], [18, 22, 2.0], [18, 23, 2.0], [18, 24, 9.0], [18, 28, 1.0], [18, 29, 1.0], [18, 30, 2.0], [18, 34, 2.0], [19, 1, 2.0], [19, 2, 2.0], [19, 3, 1.0], [19, 5, 1.0], [19, 6, 4.0], [19, 27, 1.0], [19, 36, 1.0], [20, 12, 2.0], [20, 16, 1.0], [20, 17, 1.0], [20, 18, 4.0], [20, 19, 2.0], [20, 20, 2.0], [20, 21, 3.0], [20, 24, 1.0], [21, 0, 1.0], [21, 5, 1.0], [21, 8, 2.0], [21, 9, 3.0], [22, 12, 2.0], [22, 14, 1.0], [22, 15, 2.0], [22, 16, 3.0], [22, 22, 2.0], [22, 23, 2.0], [23, 25, 1.0], [23, 26, 1.0], [23, 27, 1.0], [23, 29, 1.0], [23, 30, 2.0], [23, 32, 1.0], [23, 35, 1.0], [23, 36, 3.0], [24, 4, 1.0], [24, 6, 1.0], [24, 8, 1.0], [24, 12, 9.0], [24, 13, 3.0], [24, 14, 2.0], [24, 15, 3.0], [24, 16, 18.0], [24, 17, 3.0], [24, 18, 3.0], [24, 19, 2.0], [24, 20, 3.0], [24, 21, 5.0], [24, 22, 1.0], [24, 23, 1.0], [24, 24, 1.0], [25, 2, 1.0], [25, 5, 1.0], [25, 6, 1.0], [25, 7, 1.0], [25, 8, 1.0], [25, 9, 1.0], [25, 10, 1.0], [25, 12, 4.0], [25, 13, 3.0], [25, 14, 1.0], [25, 15, 1.0], [25, 16, 1.0], [25, 18, 2.0], [25, 20, 6.0], [25, 21, 3.0], [25, 22, 3.0], [25, 23, 1.0], [25, 24, 3.0], [25, 28, 1.0], [25, 29, 1.0], [25, 30, 1.0], [25, 32, 1.0], [25, 34, 1.0], [25, 37, 1.0], [26, 12, 2.0], [26, 13, 3.0], [26, 14, 4.0], [26, 15, 1.0], [26, 18, 3.0], [26, 19, 1.0], [26, 20, 1.0], [26, 21, 4.0], [26, 22, 1.0], [26, 23, 4.0], [26, 24, 1.0], [26, 28, 1.0], [27, 12, 2.0], [27, 14, 1.0], [27, 15, 1.0], [27, 17, 2.0], [27, 18, 2.0], [27, 21, 1.0], [27, 22, 1.0], [28, 0, 150.0], [28, 1, 149.0], [28, 2, 154.0], [28, 3, 58.0], [28, 4, 156.0], [28, 5, 137.0], [28, 6, 153.0], [28, 7, 163.0], [28, 8, 151.0], [28, 9, 157.0], [28, 10, 176.0], [28, 11, 155.0], [28, 12, 16.0], [28, 13, 17.0], [28, 14, 9.0], [28, 15, 38.0], [28, 16, 6.0], [28, 17, 8.0], [28, 18, 3.0], [28, 19, 14.0], [28, 20, 54.0], [28, 21, 37.0], [28, 22, 52.0], [28, 23, 17.0], [28, 24, 29.0], [28, 25, 137.0], [28, 26, 153.0], [28, 27, 140.0], [28, 28, 110.0], [28, 29, 132.0], [28, 30, 97.0], [28, 31, 128.0], [28, 32, 161.0], [28, 33, 177.0], [28, 34, 143.0], [28, 35, 141.0], [28, 36, 126.0], [28, 37, 123.0], [29, 12, 1.0], [29, 13, 1.0], [29, 16, 1.0], [29, 17, 1.0], [29, 18, 1.0], [29, 22, 1.0], [30, 3, 2.0], [30, 6, 1.0], [30, 12, 5.0], [30, 13, 8.0], [30, 14, 5.0], [30, 15, 3.0], [30, 16, 1.0], [30, 18, 4.0], [30, 20, 8.0], [30, 21, 9.0], [30, 22, 5.0], [30, 23, 2.0], [30, 24, 3.0], [30, 26, 1.0], [30, 27, 2.0], [30, 32, 1.0], [30, 37, 1.0], [31, 2, 1.0], [31, 3, 1.0], [31, 5, 2.0], [31, 6, 1.0], [31, 9, 3.0], [31, 10, 1.0], [31, 12, 12.0], [31, 13, 2.0], [31, 14, 2.0], [31, 15, 2.0], [31, 16, 1.0], [31, 18, 4.0], [31, 20, 5.0], [31, 21, 6.0], [31, 22, 2.0], [31, 23, 4.0], [31, 24, 4.0], [31, 25, 14.0], [31, 27, 3.0], [31, 28, 23.0], [31, 29, 18.0], [31, 30, 12.0], [31, 31, 11.0], [31, 32, 4.0], [31, 33, 1.0], [31, 34, 4.0], [31, 35, 1.0], [31, 36, 1.0], [31, 37, 4.0], [32, 5, 1.0], [32, 8, 2.0], [32, 9, 1.0], [32, 12, 5.0], [32, 14, 8.0], [32, 15, 3.0], [32, 16, 5.0], [32, 17, 43.0], [32, 18, 13.0], [32, 19, 30.0], [32, 20, 1.0], [32, 21, 1.0], [32, 22, 2.0], [32, 23, 3.0], [32, 24, 20.0], [32, 36, 1.0], [33, 12, 2.0], [33, 15, 1.0], [33, 16, 2.0], [33, 20, 2.0], [33, 21, 1.0], [33, 23, 1.0], [34, 3, 1.0], [34, 12, 2.0], [34, 13, 1.0], [34, 15, 3.0], [34, 18, 1.0], [34, 19, 1.0], [34, 20, 2.0], [34, 23, 2.0], [35, 26, 1.0], [35, 27, 1.0], [35, 28, 1.0], [35, 29, 2.0], [35, 30, 1.0], [36, 2, 1.0], [36, 3, 1.0], [36, 4, 1.0], [37, 12, 2.0], [37, 14, 1.0], [37, 18, 1.0], [37, 20, 1.0], [37, 21, 1.0], [37, 23, 1.0], [38, 1, 1.0], [38, 2, 2.0], [38, 12, 2.0], [38, 13, 8.0], [38, 14, 34.0], [38, 15, 2.0], [38, 16, 2.0], [38, 17, 2.0], [38, 18, 6.0], [38, 19, 35.0], [38, 20, 1.0], [38, 21, 4.0], [38, 22, 4.0], [38, 23, 3.0], [38, 24, 3.0], [38, 25, 1.0], [38, 31, 1.0], [39, 12, 2.0], [39, 13, 4.0], [39, 14, 1.0], [39, 17, 7.0], [39, 18, 6.0], [39, 19, 3.0], [39, 20, 2.0], [39, 23, 1.0], [39, 24, 2.0], [40, 14, 1.0], [40, 16, 3.0], [40, 17, 2.0], [40, 18, 1.0], [40, 19, 3.0], [40, 23, 2.0], [41, 2, 1.0], [41, 4, 10.0], [41, 5, 1.0], [41, 9, 2.0], [41, 11, 2.0], [41, 12, 5.0], [41, 13, 3.0], [41, 14, 1.0], [41, 15, 2.0], [41, 18, 2.0], [41, 19, 1.0], [41, 20, 4.0], [41, 21, 2.0], [41, 22, 5.0], [41, 23, 14.0], [41, 24, 7.0], [41, 25, 1.0], [41, 26, 1.0], [41, 27, 2.0], [41, 28, 4.0], [41, 30, 14.0], [41, 32, 2.0], [41, 33, 2.0], [41, 35, 2.0], [41, 36, 2.0], [41, 37, 5.0], [42, 2, 2.0], [42, 5, 1.0], [42, 9, 2.0], [42, 14, 6.0], [42, 15, 2.0], [42, 17, 9.0], [42, 19, 3.0], [42, 23, 2.0], [42, 24, 7.0], [42, 25, 3.0], [42, 26, 1.0], [42, 28, 2.0], [42, 30, 1.0], [42, 32, 2.0], [42, 33, 1.0], [42, 34, 3.0], [42, 35, 1.0], [42, 36, 1.0], [43, 12, 1.0], [43, 17, 4.0], [43, 18, 1.0], [43, 19, 2.0], [43, 23, 1.0], [43, 24, 3.0], [43, 36, 1.0], [44, 0, 2.0], [44, 1, 3.0], [44, 2, 10.0], [44, 3, 96.0], [44, 4, 2.0], [44, 5, 4.0], [44, 6, 7.0], [44, 7, 8.0], [44, 8, 2.0], [44, 9, 5.0], [44, 10, 1.0], [44, 12, 18.0], [44, 13, 13.0], [44, 14, 12.0], [44, 15, 39.0], [44, 16, 13.0], [44, 17, 4.0], [44, 18, 16.0], [44, 19, 2.0], [44, 20, 32.0], [44, 21, 21.0], [44, 22, 32.0], [44, 23, 13.0], [44, 24, 8.0], [44, 25, 1.0], [44, 27, 2.0], [44, 29, 1.0], [44, 30, 5.0], [44, 32, 3.0], [44, 34, 3.0], [44, 35, 3.0], [44, 36, 3.0], [45, 0, 1.0], [45, 2, 1.0], [45, 3, 1.0], [45, 6, 1.0], [45, 7, 1.0], [45, 9, 2.0], [45, 10, 1.0], [45, 12, 1.0], [45, 13, 1.0], [45, 15, 3.0], [45, 18, 1.0], [45, 22, 2.0], [45, 24, 2.0], [45, 25, 1.0], [45, 30, 2.0], [45, 35, 1.0], [46, 4, 1.0], [46, 12, 2.0], [46, 13, 2.0], [46, 14, 3.0], [46, 15, 3.0], [46, 16, 3.0], [46, 17, 19.0], [46, 18, 3.0], [46, 19, 9.0], [46, 21, 1.0], [46, 22, 2.0], [46, 23, 5.0], [46, 24, 4.0], [46, 26, 1.0], [46, 31, 3.0], [47, 13, 1.0], [47, 14, 1.0], [47, 16, 7.0], [47, 17, 1.0], [47, 18, 1.0], [47, 21, 2.0], [48, 12, 3.0], [48, 14, 3.0], [48, 16, 3.0], [48, 19, 1.0], [48, 21, 1.0], [48, 23, 1.0], [49, 25, 1.0], [49, 26, 1.0], [49, 28, 2.0], [49, 29, 1.0], [49, 30, 2.0], [49, 35, 1.0]], "columns": [{"id": "S1RingL", "metadata": null}, {"id": "S1keyM", "metadata": null}, {"id": "S1keySpace", "metadata": null}, {"id": "S1IndexL", "metadata": null}, {"id": "S1keyK", "metadata": null}, {"id": "S1ThumbR", "metadata": null}, {"id": "S1keyV", "metadata": null}, {"id": "S1IndexR", "metadata": null}, {"id": "S1keyA", "metadata": null}, {"id": "S1RingR", "metadata": null}, {"id": "S1MiddleR", "metadata": null}, {"id": "S1keyD", "metadata": null}, {"id": "S2keySpace", "metadata": null}, {"id": "S2keyJ", "metadata": null}, {"id": "S2keyLeftShift", "metadata": null}, {"id": "S2keyN", "metadata": null}, {"id": "S2keyZ", "metadata": null}, {"id": "S2IndexL", "metadata": null}, {"id": "S2keyA", "metadata": null}, {"id": "S2PinkyL", "metadata": null}, {"id": "S2keyK", "metadata": null}, {"id": "S2keyRightShift", "metadata": null}, {"id": "S2keyM", "metadata": null}, {"id": "S2keyI", "metadata": null}, {"id": "S2PinkyR", "metadata": null}, {"id": "S3keySpace", "metadata": null}, {"id": "S3keyEnter", "metadata": null}, {"id": "S3keyS", "metadata": null}, {"id": "S3IndexR", "metadata": null}, {"id": "S3ThumbR", "metadata": null}, {"id": "S3MiddleR", "metadata": null}, {"id": "S3keyY", "metadata": null}, {"id": "S3ThumbL", "metadata": null}, {"id": "S3keyF", "metadata": null}, {"id": "S3IndexL", "metadata": null}, {"id": "S3keyW", "metadata": null}, {"id": "S3keyQ", "metadata": null}, {"id": "S3keyL", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2595", "matrix_type": "sparse", "shape": [50, 38], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-22T01:46:08.091846", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""


test_map=\
'#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tIndividual\tDescription\rS1RingL\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1keyM\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1keySpace\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1IndexL\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1keyK\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1ThumbR\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1keyV\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1IndexR\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1keyA\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1RingR\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1MiddleR\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS1keyD\t\t\tS1\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2keySpace\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2keyJ\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2keyLeftShift\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2keyN\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2keyZ\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2IndexL\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2keyA\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2PinkyL\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2keyK\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2keyRightShift\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2keyM\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2keyI\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS2PinkyR\t\t\tS2\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3keySpace\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3keyEnter\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3keyS\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3IndexR\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3ThumbR\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3MiddleR\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3keyY\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3ThumbL\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3keyF\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3IndexL\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3keyW\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3keyQ\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"\rS3keyL\t\t\tS3\t"QIIME test code mapping file, similar to that in supervised_learning.py"'
# run unit tests if run from command-line
if __name__ == '__main__':
    main()
