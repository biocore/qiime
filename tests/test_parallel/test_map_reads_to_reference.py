#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from glob import glob
from shutil import rmtree
from os.path import exists, join
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import create_dir, remove_files
from biom.parse import parse_biom_table
from qiime.test import initiate_timeout, disable_timeout
from qiime.util import get_qiime_temp_dir, get_tmp_filename
from qiime.parse import parse_otu_map
from qiime.parallel.map_reads_to_reference import \
    (ParallelDatabaseMapperUsearch, ParallelDatabaseMapperBlat,
     ParallelDatabaseMapperBwaShort)


class ParallelDatabaseMapperTests(TestCase):

    def setUp(self):
        """ """
        self.files_to_remove = []
        self.dirs_to_remove = []

        tmp_dir = get_qiime_temp_dir()
        self.test_out = get_tmp_filename(tmp_dir=tmp_dir,
                                         prefix='qiime_parallel_tests_',
                                         suffix='',
                                         result_constructor=str)
        self.dirs_to_remove.append(self.test_out)
        create_dir(self.test_out)

        self.refseqs1_fp = get_tmp_filename(tmp_dir=self.test_out,
                                            prefix='qiime_refseqs',
                                            suffix='.fasta')
        refseqs1_f = open(self.refseqs1_fp, 'w')
        refseqs1_f.write(refseqs1)
        refseqs1_f.close()
        self.files_to_remove.append(self.refseqs1_fp)

        self.refseqs2_fp = get_tmp_filename(tmp_dir=self.test_out,
                                            prefix='qiime_refseqs',
                                            suffix='.fasta')
        refseqs2_f = open(self.refseqs2_fp, 'w')
        refseqs2_f.write(refseqs2)
        refseqs2_f.close()
        self.files_to_remove.append(self.refseqs2_fp)

        self.inseqs1_fp = get_tmp_filename(tmp_dir=self.test_out,
                                           prefix='qiime_inseqs',
                                           suffix='.fasta')
        inseqs1_f = open(self.inseqs1_fp, 'w')
        inseqs1_f.write(inseqs1)
        inseqs1_f.close()
        self.files_to_remove.append(self.inseqs1_fp)

        self.inseqs2_fp = get_tmp_filename(tmp_dir=self.test_out,
                                           prefix='qiime_inseqs',
                                           suffix='.fasta')
        inseqs2_f = open(self.inseqs2_fp, 'w')
        inseqs2_f.write(inseqs2)
        inseqs2_f.close()
        self.files_to_remove.append(self.inseqs2_fp)

        initiate_timeout(60)

    def tearDown(self):
        """ """
        disable_timeout()
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)


class ParallelDatabaseMapperUsearchTests(ParallelDatabaseMapperTests):

    def test_parallel_database_mapper_usearch(self):
        """ parallel_database_mapper_usearch functions as expected """

        params = {'refseqs_fp': self.refseqs1_fp,
                  'min_percent_id': 0.97,
                  'evalue': 1e-10,
                  'max_accepts': 1,
                  'max_rejects': 32,
                  'queryalnfract': 0.35,
                  'targetalnfract': 0.0,
                  'observation_metadata_fp': None
                  }

        app = ParallelDatabaseMapperUsearch()
        r = app(self.inseqs1_fp,
                self.test_out,
                params,
                job_prefix='PTEST',
                poll_directly=True,
                suppress_submit_jobs=False)
        observation_map_fp = glob(
            join(self.test_out, 'observation_map.txt'))[0]
        omap = parse_otu_map(open(observation_map_fp, 'U'))
        self.assertEqual(len(omap[0]), 3)
        self.assertEqualItems(
            omap[1],
            ['eco:b0015',
             'eco:b0122',
             'eco:b0015:duplicate'])
        self.assertEqualItems(omap[2], ['eco:b0015-pr', 'eco:b0122-pr'])


class ParallelDatabaseMapperBlatTests(ParallelDatabaseMapperTests):

    def test_parallel_database_mapper_blat(self):
        """ parallel_database_mapper_blat functions as expected """

        params = {'refseqs_fp': self.refseqs1_fp,
                  'min_percent_id': 0.97,
                  'evalue': 1e-10,
                  'max_accepts': 1,
                  'max_rejects': 32,
                  'queryalnfract': 0.35,
                  'targetalnfract': 0.0,
                  'observation_metadata_fp': None
                  }

        app = ParallelDatabaseMapperBlat()
        r = app(self.inseqs1_fp,
                self.test_out,
                params,
                job_prefix='PTEST',
                poll_directly=True,
                suppress_submit_jobs=False)
        observation_map_fp = glob(
            join(self.test_out, 'observation_map.txt'))[0]
        omap = parse_otu_map(open(observation_map_fp, 'U'))
        self.assertEqual(len(omap[0]), 3)
        self.assertEqualItems(
            omap[1],
            ['eco:b0015',
             'eco:b0122',
             'eco:b0015:duplicate'])
        self.assertEqualItems(omap[2], ['eco:b0015-pr', 'eco:b0122-pr'])


class ParallelDatabaseMapperBwaShortTests(ParallelDatabaseMapperTests):

    def test_bwa_short_database_mapper(self):
        """bwa_short_database_mapper functions as expected """
        params = {'refseqs_fp': self.refseqs2_fp,
                  'max_diff': None,
                  'observation_metadata_fp': None}
        app = ParallelDatabaseMapperBwaShort()
        r = app(self.inseqs2_fp,
                self.test_out,
                params,
                poll_directly=True,
                suppress_submit_jobs=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')
        table = parse_biom_table(open(observation_table_fp, 'U'))
        self.assertEqualItems(table.SampleIds, ['s2', 's1'])
        self.assertEqualItems(
            table.ObservationIds,
            ['r1',
             'r2',
             'r3',
             'r4',
             'r5'])
        self.assertEqual(table.sum(), 6)

    def test_bwa_short_database_mapper_alt_params(self):
        """bwa_short_database_mapper functions as expected """
        params = {'refseqs_fp': self.refseqs2_fp,
                  'max_diff': 1,
                  'observation_metadata_fp': None}
        app = ParallelDatabaseMapperBwaShort()
        r = app(self.inseqs2_fp,
                self.test_out,
                params,
                poll_directly=True,
                suppress_submit_jobs=False)
        observation_map_fp = join(self.test_out, 'observation_map.txt')
        self.assertTrue(exists(observation_map_fp))
        observation_table_fp = join(self.test_out, 'observation_table.biom')
        table = parse_biom_table(open(observation_table_fp, 'U'))
        self.assertEqualItems(table.SampleIds, ['s2', 's1'])
        self.assertEqualItems(table.ObservationIds, ['r2', 'r3', 'r4', 'r5'])
        self.assertEqual(table.sum(), 5)

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

inseqs1 = """>eco:b0001 thrL; thr operon leader peptide; K08278 thr operon leader peptide (N)
atgaaacgcattagcaccaccattaccaccaccatcaccattaccacaggtaacggtgcg
ggctga
>eco:b0015 dnaJ; chaperone Hsp40, co-chaperone with DnaK; K03686 molecular chaperone DnaJ (N)
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
>eco:b0122
atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaataa
>eco:b0015:duplicate
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
"""

inseqs2 = """>s1_1
atgttacgcattagcaccaccattaccaccaccatcaccattaccacaggtaacggtgcg
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
