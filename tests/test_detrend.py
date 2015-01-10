#!/usr/bin/env python

"""Tests of code for detrending"""

__author__ = "Dan Knights"
__copyright__ = "Copyright 2012, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.9.0-rc2"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"


from os import remove, system, close
from shutil import rmtree
from os.path import join, exists
from tempfile import NamedTemporaryFile, mkdtemp, mkstemp

from unittest import TestCase, main
from burrito.util import ApplicationError

from skbio.util import remove_files
from qiime.detrend import detrend_pcoa
from numpy import array


def is_float(input_string):
    """True if string can be cast as a float"""
    try:
        float(input_string)
        return True
    except ValueError:
        return False


class DetrendTests(TestCase):

    """Tests of the RSupervisedLearner class"""

    def setUp(self):

        # Temporary input file
        fd, self.tmp_pc_fp = mkstemp(prefix='R_test_pcoa',
                                    suffix='.txt')
        close(fd)
        seq_file = open(self.tmp_pc_fp, 'w')
        seq_file.write(test_pc)
        seq_file.close()

        fd, self.tmp_map_fp = mkstemp(prefix='R_test_map_',
                                     suffix='.txt')
        close(fd)
        map_file = open(self.tmp_map_fp, 'w')
        map_file.write(test_map)
        map_file.close()

        self.files_to_remove = \
            [self.tmp_pc_fp, self.tmp_map_fp]

        # Prep input files in R format
        self.output_dir = mkdtemp()
        self.dirs_to_remove = [self.output_dir]

        # get detrending results
#         mkdir(join(self.output_dir, 'random_forest'))
#         self.results = RSupervisedLearner()(
#             self.tmp_otu_fp, self.tmp_map_fp,'Individual',
#             ntree=100, errortype='oob',
#             output_dir=self.output_dir)

    def tearDown(self):
        return
        remove_files(set(self.files_to_remove))
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_detrend_no_gradient(self):
        """Ensure that detrending runs and produces expected files
           with no associated gradient.

        """
        results = detrend_pcoa(input_fp=self.tmp_pc_fp,
                               map_fp=None, gradient_variable=None,
                               suppress_prerotate=False, output_dir=self.output_dir,
                               HALT_EXEC=False)
        self.assertEqual(results['summary'], None)
        coords = results['coords']
        lines = coords.readlines()

        # ensure one line per sample in detrended pcoa
        self.assertEqual(len(lines), len(test_pc.split('\n')) - 4)
        # ensure three columns tab delimited
        self.assertEqual(len(lines[0].split('\t')), 3)

    def test_detrend_gradient(self):
        """Ensure that detrending runs and produces expected files
           with associated gradient.

        """
        results = detrend_pcoa(input_fp=self.tmp_pc_fp,
                               map_fp=self.tmp_map_fp, gradient_variable='Gradient',
                               suppress_prerotate=False, output_dir=self.output_dir,
                               HALT_EXEC=False)

        # check formatting of summary file
        lines = results['summary'].readlines()
        self.assertEqual(len(lines), 4)

        # check formatting of coords file
        lines = results['coords'].readlines()
        # ensure one line per sample in detrended pcoa
        self.assertEqual(len(lines), len(test_pc.split('\n')) - 4)
        # ensure three columns tab delimited
        self.assertEqual(len(lines[0].split('\t')), 3)

        # ensure that plot pdf is at least present
        self.assertEqual(str(type(results['plot'])), "<type 'file'>")

    def test_detrend_gradient_no_prerotate(self):
        """Ensure that detrending runs and produces expected files
           with associated gradient.

        """
        results = detrend_pcoa(input_fp=self.tmp_pc_fp,
                               map_fp=self.tmp_map_fp, gradient_variable='Gradient',
                               suppress_prerotate=True, output_dir=self.output_dir,
                               HALT_EXEC=False)

        # check formatting of summary file
        lines = results['summary'].readlines()
        self.assertEqual(len(lines), 4)

        # check formatting of coords file
        lines = results['coords'].readlines()
        # ensure one line per sample in detrended pcoa
        self.assertEqual(len(lines), len(test_pc.split('\n')) - 4)
        # ensure three columns tab delimited
        self.assertEqual(len(lines[0].split('\t')), 3)

        # ensure that plot pdf is at least present
        self.assertEqual(str(type(results['plot'])), "<type 'file'>")


test_pc = \
    """pc vector number	1	2
outsidemouthT1.U1.536668	-0.259989511	0.031981421
outsidemouthT1.U2.536736	-0.10770442	0.00760908
outsidemouthT2.U1.536267	0.007775311	0.085673368
outsidemouthT2.U2.536494	0.096648247	0.11165307
outsidemouthT3.U1.536512	0.060330843	0.097689673
outsidemouthT3.U2.536209	0.086096188	0.120917046
outsidemouthT4.U1.536532	0.104765697	-0.005655307
outsidemouthT4.U2.536258	0.154786488	0.146644973
outsidemouthT5.U1.536645	0.111274125	-0.016418682
outsidemouthT5.U2.536511	0.09092745	-0.03218781
outsidemouthT6.U1.536444	0.070454432	-0.212724803
outsidemouthT6.U2.536486	0.067498055	-0.241014047
outsidemouthT7.U1.536739	0.030952249	-0.260790565
outsidemouthT7.U2.536211	0.031283845	-0.262058153

eigvals	2.231000323	1.471086217
% variation explained	22.05856065	14.54506492
"""

test_map = \
    """#SampleID	Gradient
outsidemouthT1.U1.536668	10
outsidemouthT1.U2.536736	20
outsidemouthT2.U1.536267	15
outsidemouthT2.U2.536494	25
outsidemouthT3.U1.536512	10
outsidemouthT3.U2.536209	10
outsidemouthT4.U1.536532	30
outsidemouthT4.U2.536258	35
outsidemouthT5.U1.536645	1000
outsidemouthT5.U2.536511	18
outsidemouthT6.U1.536444	17
outsidemouthT6.U2.536486	16
outsidemouthT7.U1.536739	12
outsidemouthT7.U2.536211	11
"""

# run unit tests if run from command-line
if __name__ == '__main__':
    main()
