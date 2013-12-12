#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Michael Dwan", "Logan Knecht",
               "Damien Coy", "Levi McCracken"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the compare_categories.py module."""

from os.path import exists, join
from shutil import rmtree
from tempfile import mkdtemp, NamedTemporaryFile
from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from qiime.compare_categories import compare_categories
from qiime.util import get_qiime_temp_dir

class CompareCategoriesTests(TestCase):
    """Tests for the compare_categories.py module."""

    def setUp(self):
        """Defines data that will be used by the tests."""
        self.files_to_remove = []
        self.dirs_to_remove = []

        # Create temp directory to hold input and output.
        self.test_dir = mkdtemp(dir=get_qiime_temp_dir(),
                                prefix='qiime_compare_categories_tests_')
        self.dirs_to_remove.append(self.test_dir)

        # Create input files under our temp dir.
        self.dm1_fp = join(self.test_dir, 'dm1.txt')
        dm1_f = open(self.dm1_fp, 'w')
        dm1_f.write(dm1_str)
        dm1_f.close()
        self.files_to_remove.append(self.dm1_fp)

        self.dm2_fp = join(self.test_dir, 'dm2.txt')
        dm2_f = open(self.dm2_fp, 'w')
        dm2_f.write(dm2_str)
        dm2_f.close()
        self.files_to_remove.append(self.dm2_fp)

        self.invalid_dm_fp = join(self.test_dir, 'invalid_dm.txt')
        invalid_dm_f = open(self.invalid_dm_fp, 'w')
        invalid_dm_f.write(invalid_dm_str)
        invalid_dm_f.close()
        self.files_to_remove.append(self.invalid_dm_fp)

        self.map1_fp = join(self.test_dir, 'map1.txt')
        map1_f = open(self.map1_fp, 'w')
        map1_f.write(map1_str)
        map1_f.close()
        self.files_to_remove.append(self.map1_fp)

        self.map2_fp = join(self.test_dir, 'map2.txt')
        map2_f = open(self.map2_fp, 'w')
        map2_f.write(map2_str)
        map2_f.close()
        self.files_to_remove.append(self.map2_fp)

        self.missing_sample_map_fp = join(self.test_dir,
                                          'missing_sample_map_fp.txt')
        missing_sample_map_f = open(self.missing_sample_map_fp, 'w')
        missing_sample_map_f.write('\n'.join(map1_str.split('\n')[:-3]))
        missing_sample_map_f.close()
        self.files_to_remove.append(self.missing_sample_map_fp)

        self.cat_methods = ['adonis', 'anosim', 'mrpp', 'permanova',
                            'permdisp', 'dbrda']
        self.num_methods = ['best', 'morans_i']
        self.cat_categories = ['Treatment']
        self.num_categories = ['DOB']
        self.num_perms = 42

    def tearDown(self):
        """Removes temporary directories and files."""
        remove_files(self.files_to_remove)

        # Remove directories last, so we don't get errors trying to remove
        # files which may be in the directories.
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_compare_categories_categorical_variables(self):
        """Test compare_categories() on categorical methods/input."""
        for method in self.cat_methods:
            compare_categories(self.dm1_fp, self.map1_fp, method,
                    self.cat_categories, self.num_perms, self.test_dir)
            results_fp = join(self.test_dir, '%s_results.txt' % method)
            self.files_to_remove.append(results_fp)
            results_f = open(results_fp, 'U')
            results = results_f.readlines()
            results_f.close()

            # Make sure the files aren't empty.
            self.assertTrue(len(results) > 0)

    def test_compare_categories_numeric_variables(self):
        """Test compare_categories() on numeric methods/input."""
        for method in self.num_methods:
            compare_categories(self.dm1_fp, self.map1_fp, method,
                    self.num_categories, self.num_perms, self.test_dir)
            results_fp = join(self.test_dir, '%s_results.txt' % method)
            self.files_to_remove.append(results_fp)
            results_f = open(results_fp, 'U')
            results = results_f.readlines()
            results_f.close()
            self.assertTrue(len(results) > 0)

    def test_compare_categories_morans_i_zeros(self):
        """Test Moran's I on distance matrix with non-diagonal zeros."""
        method = 'morans_i'
        compare_categories(self.dm2_fp, self.map2_fp, method, ['NumCat'], 999,
                           self.test_dir)
        results_fp = join(self.test_dir, '%s_results.txt' % method)
        self.files_to_remove.append(results_fp)
        results_f = open(results_fp, 'U')
        results = results_f.readlines()
        results_f.close()
        self.assertTrue(len(results) > 0)

    def test_compare_categories_invalid_input(self):
        """Test compare_categories() on invalid input that should error out."""
        # Non-numeric categories with BEST and Moran's I.
        for method in self.num_methods:
            self.assertRaises(TypeError, compare_categories, self.dm1_fp,
                    self.map1_fp, method, self.cat_categories, self.num_perms,
                    self.test_dir)

        # SampleID with all methods.
        for method in self.num_methods + self.cat_methods:
            self.assertRaises(ValueError, compare_categories, self.dm1_fp,
                    self.map1_fp, method, ['SampleID'], self.num_perms,
                    self.test_dir)

        # Single category passed as a string instead of a list of string(s).
        for method in self.num_methods + self.cat_methods:
            self.assertRaises(TypeError, compare_categories, self.dm1_fp,
                    self.map1_fp, method, 'SampleID', self.num_perms,
                    self.test_dir)

        # Non-symmetric/hollow distance matrix.
        for method in self.num_methods:
            self.assertRaises(ValueError, compare_categories,
                    self.invalid_dm_fp, self.map1_fp, method,
                    self.num_categories, self.num_perms, self.test_dir)
        for method in self.cat_methods:
            self.assertRaises(ValueError, compare_categories,
                    self.invalid_dm_fp, self.map1_fp, method,
                    self.cat_categories, self.num_perms, self.test_dir)

        # Nonexistent category.
        for method in self.num_methods + self.cat_methods:
            self.assertRaises(ValueError, compare_categories, self.dm1_fp,
                    self.map1_fp, method, ['bar'], self.num_perms,
                    self.test_dir)

        # Unique category values only.
        for method in self.cat_methods:
            self.assertRaises(ValueError, compare_categories, self.dm1_fp,
                    self.map1_fp, method, ['Unique'], self.num_perms,
                    self.test_dir)

        # Only a single category value.
        for method in self.cat_methods + self.num_methods:
            if method == 'best':
                # BEST is okay with this type of category.
                compare_categories(self.dm1_fp,
                        self.map1_fp, method, ['Single'], self.num_perms,
                        self.test_dir)
                results_fp = join(self.test_dir, '%s_results.txt' % method)
                self.files_to_remove.append(results_fp)
                results_f = open(results_fp, 'U')
                results = results_f.readlines()
                results_f.close()
                self.assertTrue(len(results) > 0)
            else:
                self.assertRaises(ValueError, compare_categories, self.dm1_fp,
                        self.map1_fp, method, ['Single'], self.num_perms,
                        self.test_dir)

        # Bad number of permutations.
        for method in self.cat_methods:
            self.assertRaises(ValueError, compare_categories, self.dm1_fp,
                    self.map1_fp, method, self.cat_categories, -42,
                    self.test_dir)

        # Unrecognized method.
        self.assertRaises(ValueError, compare_categories, self.dm1_fp,
                          self.map1_fp, 'foo', self.cat_categories,
                          self.num_perms, self.test_dir)

        # Samples in dm not found in map.
        for method in self.num_methods:
            with self.assertRaises(ValueError):
                compare_categories(self.dm1_fp, self.missing_sample_map_fp,
                                   method, self.num_categories, self.num_perms,
                                   self.test_dir)
        for method in self.cat_methods:
            with self.assertRaises(ValueError):
                compare_categories(self.dm1_fp, self.missing_sample_map_fp,
                                   method, self.cat_categories, self.num_perms,
                                   self.test_dir)


dm1_str = """\tPC.354\tPC.355\tPC.356\tPC.481\tPC.593\tPC.607\tPC.634\tPC.635\tPC.636
PC.354\t0.0\t0.595483768391\t0.618074717633\t0.582763100909\t0.566949022108\t0.714717232268\t0.772001731764\t0.690237118413\t0.740681707488
PC.355\t0.595483768391\t0.0\t0.581427669668\t0.613726772383\t0.65945132763\t0.745176523638\t0.733836123821\t0.720305073505\t0.680785600439
PC.356\t0.618074717633\t0.581427669668\t0.0\t0.672149021573\t0.699416863323\t0.71405573754\t0.759178215168\t0.689701276341\t0.725100672826
PC.481\t0.582763100909\t0.613726772383\t0.672149021573\t0.0\t0.64756120797\t0.666018240373\t0.66532968784\t0.650464714994\t0.632524644216
PC.593\t0.566949022108\t0.65945132763\t0.699416863323\t0.64756120797\t0.0\t0.703720200713\t0.748240937349\t0.73416971958\t0.727154987937
PC.607\t0.714717232268\t0.745176523638\t0.71405573754\t0.666018240373\t0.703720200713\t0.0\t0.707316869557\t0.636288883818\t0.699880573956
PC.634\t0.772001731764\t0.733836123821\t0.759178215168\t0.66532968784\t0.748240937349\t0.707316869557\t0.0\t0.565875193399\t0.560605525642
PC.635\t0.690237118413\t0.720305073505\t0.689701276341\t0.650464714994\t0.73416971958\t0.636288883818\t0.565875193399\t0.0\t0.575788039321
PC.636\t0.740681707488\t0.680785600439\t0.725100672826\t0.632524644216\t0.727154987937\t0.699880573956\t0.560605525642\t0.575788039321\t0.0
"""

# Some zeros in the matrix (not just on the diagonal).
dm2_str = """\tS1\tS2\tS3
S1\t0.0\t0.0\t0.5
S2\t0.0\t0.0\t0.1
S3\t0.5\t0.1\t0.0
"""

# Non-symmetric. :(
invalid_dm_str = """\tPC.354\tPC.355\tPC.356\tPC.481\tPC.593\tPC.607\tPC.634\tPC.635\tPC.636
PC.354\t0.0\t0.599483768391\t0.618074717633\t0.582763100909\t0.566949022108\t0.714717232268\t0.772001731764\t0.690237118413\t0.740681707488
PC.355\t0.595483768391\t0.0\t0.581427669668\t0.613726772383\t0.65945132763\t0.745176523638\t0.733836123821\t0.720305073505\t0.680785600439
PC.356\t0.618074717633\t0.581427669668\t0.0\t0.672149021573\t0.699416863323\t0.71405573754\t0.759178215168\t0.689701276341\t0.725100672826
PC.481\t0.582763100909\t0.613726772383\t0.672149021573\t0.0\t0.64756120797\t0.666018240373\t0.66532968784\t0.650464714994\t0.632524644216
PC.593\t0.566949022108\t0.65945132763\t0.699416863323\t0.64756120797\t0.0\t0.703720200713\t0.748240937349\t0.73416971958\t0.727154987937
PC.607\t0.714717232268\t0.745176523638\t0.71405573754\t0.666018240373\t0.703720200713\t0.0\t0.707316869557\t0.636288883818\t0.699880573956
PC.634\t0.772001731764\t0.733836123821\t0.759178215168\t0.66532968784\t0.748240937349\t0.707316869557\t0.0\t0.565875193399\t0.560605525642
PC.635\t0.690237118413\t0.720305073505\t0.689701276341\t0.650464714994\t0.73416971958\t0.636288883818\t0.565875193399\t0.0\t0.575788039321
PC.636\t0.740681707488\t0.680785600439\t0.725100672826\t0.632524644216\t0.727154987937\t0.699880573956\t0.560605525642\t0.575788039321\t0.0
"""

map1_str = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tDOB\tUnique\tSingle\tDescription
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354\tAGCACGAGCCTA\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\t1\t2\tControl_mouse_I.D._354
PC.355\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\t2\t2\tControl_mouse_I.D._355
PC.356\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061126\t3\t2\tControl_mouse_I.D._356
PC.481\tACCAGCGACTAG\tYATGCTGCCTCCCGTAGGAGT\tControl\t20070314\t4\t2\tControl_mouse_I.D._481
PC.593\tAGCAGCACTTGT\tYATGCTGCCTCCCGTAGGAGT\tControl\t20071210\t5\t2\tControl_mouse_I.D._593
PC.607\tAACTGTGCGTAC\tYATGCTGCCTCCCGTAGGAGT\tFast\t20071112\t6\t2\tFasting_mouse_I.D._607
PC.634\tACAGAGTCGGCT\tYATGCTGCCTCCCGTAGGAGT\tFast\t20080116\t7\t2\tFasting_mouse_I.D._634
PC.635\tACCGCAGAGTCA\tYATGCTGCCTCCCGTAGGAGT\tFast\t20080116\t8\t2\tFasting_mouse_I.D._635
PC.636\tACGGTGAGTGTC\tYATGCTGCCTCCCGTAGGAGT\tFast\t20080116\t9\t2\tFasting_mouse_I.D._636
random.sample\tACGGTGAGTGTC\tYATGCTGCCTCCCGTAGGAGT\tFast\t20080116\t9\t3\trandom.sample
"""

map2_str = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tNumCat\tDescription
S1\tAGCACGAGCCTA\tYATGCTGCCTCCCGTAGGAGT\t2\tS1
S2\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\t2\tS2
S3\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\t3\tS3
"""


if __name__ == "__main__":
    main()
