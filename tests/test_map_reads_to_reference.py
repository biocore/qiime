#!/usr/bin/env python
# File created on 13 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from os import close
from shutil import rmtree
from os.path import exists, join
from tempfile import mkstemp , mkdtemp

from skbio.util import remove_files
from unittest import TestCase, main
from numpy.testing import assert_almost_equal
from biom import load_table

from qiime.util import get_qiime_temp_dir
from qiime.map_reads_to_reference import (
    usearch_database_mapper, blat_database_mapper, bwa_short_database_mapper,
    bwa_sw_database_mapper, blat_nt_database_mapper)
from qiime.parse import parse_otu_map
from qiime.test import initiate_timeout, disable_timeout


class DatabaseAssignmentTests(TestCase):

    """
    """

    def setUp(self):
        """
        """
        tmp_dir = get_qiime_temp_dir()
        self.test_out = mkdtemp(dir=tmp_dir,
                                prefix='qiime_parallel_tests_',
                                suffix='')
        self.dirs_to_remove = [self.test_out]

        self.output_fp = join(self.test_out, 'fmap.txt')
        self.failure_fp = join(self.test_out, 'fail.txt')
        self.usearch_fp = join(self.test_out, 'out.uc')
        self.bl6_fp = join(self.test_out, 'out.bl6')
        self.log_fp = join(self.test_out, 'fmap.log')
        self.files_to_remove = [self.output_fp, self.failure_fp,
                                self.usearch_fp, self.log_fp, self.bl6_fp]

        fd, self.refseqs1_fp = mkstemp(dir=self.test_out,
                                      prefix='qiime_refseqs',
                                      suffix='.fasta')
        close(fd)
        refseqs1_f = open(self.refseqs1_fp, 'w')
        refseqs1_f.write(refseqs1)
        refseqs1_f.close()
        self.files_to_remove.append(self.refseqs1_fp)

        fd, self.refseqs2_fp = mkstemp(dir=self.test_out,
                                      prefix='qiime_refseqs',
                                      suffix='.fasta')
        close(fd)
        refseqs2_f = open(self.refseqs2_fp, 'w')
        refseqs2_f.write(refseqs2)
        refseqs2_f.close()
        self.files_to_remove.append(self.refseqs2_fp)

        fd, self.inseqs1_fp = mkstemp(dir=self.test_out,
                                     prefix='qiime_inseqs',
                                     suffix='.fasta')
        close(fd)
        inseqs1_f = open(self.inseqs1_fp, 'w')
        inseqs1_f.write(inseqs1)
        inseqs1_f.close()
        self.files_to_remove.append(self.inseqs1_fp)
        fd, self.inseqs2_fp = mkstemp(dir=self.test_out,
                                     prefix='qiime_inseqs',
                                     suffix='.fasta')
        close(fd)
        inseqs2_f = open(self.inseqs2_fp, 'w')
        inseqs2_f.write(inseqs2)
        inseqs2_f.close()
        self.files_to_remove.append(self.inseqs2_fp)
        initiate_timeout(60)

    def tearDown(self):
        """ """
        disable_timeout()
        remove_files(self.files_to_remove, error_on_missing=False)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)


class UsearchDatabaseAssignmentTests(DatabaseAssignmentTests):

    def test_usearch_database_mapper(self):
        """usearch_database_mapper functions as expected """
        usearch_database_mapper(query_fp=self.inseqs1_fp,
                                refseqs_fp=self.refseqs1_fp,
                                output_dir=self.test_out,
                                evalue=1e-10,
                                min_id=0.75,
                                queryalnfract=0.35,
                                targetalnfract=0.0,
                                maxaccepts=1,
                                maxrejects=8,
                                HALT_EXEC=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')
        table = load_table(observation_table_fp)
        self.assertItemsEqual(table.ids(), ['s2', 's1'])
        self.assertItemsEqual(
            table.ids(axis='observation'),
            ['eco:b0122-pr',
             'eco:b0015-pr'])
        self.assertEqual(table.sum(), 5)


class BlatDatabaseAssignmentTests(DatabaseAssignmentTests):

    def test_blat_database_mapper(self):
        """blat_database_mapper functions as expected """
        blat_database_mapper(query_fp=self.inseqs1_fp,
                             refseqs_fp=self.refseqs1_fp,
                             output_dir=self.test_out,
                             evalue=1e-10,
                             min_id=0.75,
                             genetic_code=11,
                             HALT_EXEC=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')
        table = load_table(observation_table_fp)
        self.assertItemsEqual(table.ids(), ['s2', 's1'])
        self.assertItemsEqual(
            table.ids(axis='observation'),
            ['eco:b0122-pr',
             'eco:b0015-pr'])
        self.assertEqual(table.sum(), 5)

    def test_blat_database_mapper_alt_params(self):
        """blat_database_mapper functions as expected w alt e-value"""
        blat_database_mapper(query_fp=self.inseqs1_fp,
                             refseqs_fp=self.refseqs1_fp,
                             output_dir=self.test_out,
                             evalue=1e-2,
                             min_id=0.75,
                             genetic_code=2,
                             HALT_EXEC=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')
        table = load_table(observation_table_fp)
        self.assertItemsEqual(table.ids(), ['s2', 's1'])
        self.assertItemsEqual(table.ids(axis='observation'),
                              ['eco:b0122-pr', 'eco:b0015-pr', 'eco:b0001-pr'])
        self.assertEqual(table.sum(), 6)


class BlatNtAssignmentTests(DatabaseAssignmentTests):

    def test_blat_nt_database_mapper(self):
        """blat_nt_database_mapper functions as expected """
        blat_nt_database_mapper(query_fp=self.inseqs2_fp,
                                refseqs_fp=self.refseqs2_fp,
                                output_dir=self.test_out,
                                evalue=1e-10,
                                min_id=0.75,
                                HALT_EXEC=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')
        table = load_table(observation_table_fp)
        self.assertItemsEqual(table.ids(), ['s2', 's1'])
        self.assertItemsEqual(
            table.ids(axis='observation'),
            ['r1',
             'r2',
             'r3',
             'r4',
             'r5'])
        self.assertEqual(table.sum(), 6)

    def test_blat_nt_database_mapper_alt_min_id(self):
        """blat_nt_database_mapper functions as expected """
        blat_nt_database_mapper(query_fp=self.inseqs2_fp,
                                refseqs_fp=self.refseqs2_fp,
                                output_dir=self.test_out,
                                evalue=1e-10,
                                min_id=1.0,
                                HALT_EXEC=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')
        table = load_table(observation_table_fp)
        self.assertItemsEqual(table.ids(), ['s2', 's1'])
        self.assertItemsEqual(table.ids(axis='observation'),
                              ['r2', 'r3', 'r4', 'r5'])
        self.assertEqual(table.sum(), 5)


class BwaShortAssignmentTests(DatabaseAssignmentTests):

    def test_bwa_short_database_mapper(self):
        """bwa_short_database_mapper functions as expected """
        bwa_short_database_mapper(query_fp=self.inseqs2_fp,
                                  refseqs_fp=self.refseqs2_fp,
                                  output_dir=self.test_out,
                                  max_diff=None,
                                  HALT_EXEC=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')
        table = load_table(observation_table_fp)
        self.assertItemsEqual(table.ids(), ['s2', 's1'])
        self.assertItemsEqual(
            table.ids(axis='observation'),
            ['r1',
             'r2',
             'r3',
             'r4',
             'r5'])
        self.assertEqual(table.sum(), 6)

    def test_bwa_short_database_mapper_alt_params(self):
        """bwa_short_database_mapper functions as expected """
        bwa_short_database_mapper(query_fp=self.inseqs2_fp,
                                  refseqs_fp=self.refseqs2_fp,
                                  output_dir=self.test_out,
                                  max_diff=1,
                                  HALT_EXEC=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')
        table = load_table(observation_table_fp)
        self.assertItemsEqual(table.ids(), ['s2', 's1'])
        self.assertItemsEqual(table.ids(axis='observation'),
                              ['r2', 'r3', 'r4', 'r5'])
        self.assertEqual(table.sum(), 5)
        # float can also be passed for max_diff
        bwa_short_database_mapper(query_fp=self.inseqs2_fp,
                                  refseqs_fp=self.refseqs2_fp,
                                  output_dir=self.test_out,
                                  max_diff=0.01,
                                  HALT_EXEC=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')


class BwaSwAssignmentTests(DatabaseAssignmentTests):

    def test_bwa_sw_database_mapper(self):
        """bwa_sw_database_mapper functions as expected """
        bwa_sw_database_mapper(query_fp=self.inseqs1_fp,
                               refseqs_fp=self.refseqs2_fp,
                               output_dir=self.test_out,
                               HALT_EXEC=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')
        table = load_table(observation_table_fp)
        self.assertItemsEqual(table.ids(), ['s2', 's1'])
        self.assertItemsEqual(
            table.ids(axis='observation'),
            ['r1',
             'r2',
             'r3',
             'r4',
             'r5'])
        self.assertEqual(table.sum(), 6)


refseqs1 = """>eco:b0001-pr
MKRISTTITTTITITTGNGAG
>eco:b0015-pr dnaJ
MAKQDYYEILGVSKTAEEREIRKAYKRLAMKYHPDRNQGDKEAEAKFKEIKEAYEVLTDS
QKRAAYDQYGHAAFEQGGMGGGGFGGGADFSDIFGDVFGDIFGGGRGRQRAARGADLRYN
MELTLEEAVRGVTKEIRIPTLEECDVCHGSGAKPGTQPQTCPTCHGSGQVQMRQGFFAVQ
QTCPHCQGRGTLIKDPCNKCHGHGRVERSKTLSVKIPAGVDTGDRIRLAGEGEAGEHGAP
AGDLYVQVQVKQHPIFEREGNNLYCEVPINFAMAALGGEIEVPTLDGRVKLKVPGETQTG
KLFRMRGKGVKSVRGGAQGDLLCRVVVETPVGLNERQKQLLQELQESFGGPTGEHNSPRS
KSFFDGVKKFFDDLTR
>eco:b0122-pr
MKTFFRTVLFGSLMAVCANSYALSESEAEDMADLTAVFVFLKNDCGYQNLPNGQIRRALV
FFAQQNQWDLSNYDTFDMKALGEDSYRDLSGIGIPVAKKCKALARDSLSLLAYVK
"""

refseqs2 = """>r1
atgaaacgcattagcaccaccattaccaccaccatcaccattaccacaggtaacggtgcg
ggctga
>r2 some comments...
atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgagaa
>r3
atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaataa
>r4
atgaagaaaattttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaatcc
>r5 some comments...
aatgactaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa
"""

inseqs1 = """>s1_1
atgaaacgcattagcaccaccattaccaccaccatcaccattaccacaggtaacggtgcg
ggctga
>s2_2 some comments...
atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgagaa
>s1_3
atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaataa
>s1_4
atgaagaaaattttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaatcc
>s1_5
atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa
>s1_6 some comments...
aatgactaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa
"""

inseqs2 = """>s1_1
atgaaacgcattagcaccaccattaccaccattatcaccattaccacaggtaacggtgcg
ggctga
>s2_2 some comments...
atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgagaa
>s1_3
atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaataa
>s1_4
atgaagaaaattttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaatcc
>s1_5
atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa
>s1_6 some comments...
aatgactaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa
"""


if __name__ == "__main__":
    main()
