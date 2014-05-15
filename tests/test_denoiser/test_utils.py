#!/usr/bin/env python

"""Tests for function in utils."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from os import remove, rmdir
from shutil import rmtree
from os.path import exists
from tempfile import mkdtemp

from unittest import TestCase, main
from numpy.testing import assert_almost_equal

from cogent import Sequence
from skbio.parse.sequences import parse_fasta
from brokit.denoiser import Flowgram, FlowgramCollection
from cogent.app.util import ApplicationNotFoundError
from cogent.util.misc import remove_files, create_dir

from qiime.util import get_qiime_project_dir
from qiime.denoiser.utils import make_stats, get_representatives,\
    squeeze_seq, wait_for_file, wait_for_cluster_ids, sort_ids,\
    sort_seqs_by_clustersize, get_denoiser_data_dir,\
    init_flowgram_file, append_to_flowgram_file, store_mapping,\
    store_clusters, invert_mapping, read_denoiser_mapping,\
    cat_sff_files, FlowgramContainerFile, FlowgramContainerArray,\
    write_checkpoint, read_checkpoint


class TestUtils(TestCase):

    def setUp(self):
        self.data = dict({"0": "ab", "1": "abababa", "2": "abab",
                          "3": "baba", "4": "ababaa", "5": "a", "6": "abababa",
                          "7": "bab", "8": "babba"})
        self.mapping = {"1": ["0", "2", "5", "6"],
                        "3": [],
                        "4": [],
                        "8": ["7"]}
        self.test_map = {'1': ('a', 'b', 'c'),
                         '2': ('d', 'e', 'f')}

        # realistic test file
        self.tiny_test = get_qiime_project_dir() +\
            "/qiime/support_files/denoiser/TestData/tiny_test.sff.txt"

        # set up test file
        open("/tmp/denoiser_utils_dummy.tmp", "w")
        self.files_to_remove = ["/tmp/denoiser_utils_dummy.tmp"]
        self.tmpdir = ""

    def tearDown(self):
        """Clean up tmp files."""
        remove_files(self.files_to_remove, False)
        if self.tmpdir:
            rmtree(self.tmpdir)

        # clean up the file from init_flowgram_file
        if (hasattr(self, "tmp_filename") and exists(self.tmp_filename)):
            remove(self.tmp_filename)

    def test_get_denoiser_data_dir(self):
        """get_denoiser_data_dir returns dir with error profiles"""

        obs = get_denoiser_data_dir()

        self.assertTrue(exists(obs))
        self.assertTrue(exists(obs + 'FLX_error_profile.dat'))

    def test_invert_mapping(self):
        """invert_prefix_map inverts a dictionary mapping."""

        actual = invert_mapping(self.test_map)
        self.assertEqual(
            {'1': '1',
             'a': '1',
             'b': '1',
             'c': '1',
             '2': '2',
             'd': '2',
             'e': '2',
             'f': '2'},
            actual)

    def test_make_stats(self):
        """make_stats produces meaningful statistics."""
        map = self.mapping
        stats = """Clustersize\t#
1:\t\t2
2:\t\t1
5:\t\t1"""

        self.assertEqual(make_stats(map), stats)

    def test_store_mapping(self):
        """store_mapping writes mapping to file."""

        expected = ["1:\t0\t2\t5\t6\n",
                    "3:\n",
                    "4:\n",
                    "8:\t7\n"]

        self.files_to_remove.append("/tmp/test_store_mapping_mapping.txt")
        store_mapping(self.mapping, "/tmp/", prefix="test_store_mapping")
        observed = list(open("/tmp/test_store_mapping_mapping.txt", "U"))
        self.assertItemsEqual(observed, expected)

    def test_store_cluster(self):
        """store_clusters stores the centroid seqs for each cluster."""

        self.tmpdir = mkdtemp(dir="./", suffix="_store_clusters/")

        self.files_to_remove.append(self.tmpdir + "singletons.fasta")
        self.files_to_remove.append(self.tmpdir + "centroids.fasta")

        # empty map results in empty files
        store_clusters({}, self.tiny_test, self.tmpdir)
        actual_centroids = list(
            parse_fasta(open(self.tmpdir + "centroids.fasta")))
        self.assertEqual(actual_centroids, [])
        actual_singletons = list(
            parse_fasta(open(self.tmpdir + "singletons.fasta")))
        self.assertEqual(actual_singletons, [])

        # non-empty map creates non-empty files, centroids sorted by size
        mapping = {'FZTHQMS01B8T1H': [],
                   'FZTHQMS01DE1KN': ['FZTHQMS01EHAJG'],
                   'FZTHQMS01EHAJG': [1, 2, 3]}  # content doesn't really matter

        centroids = [(
            'FZTHQMS01EHAJG | cluster size: 4', 'CATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGGTTTGGTGAGCCGTTACCTCACCAACTGCCTAATGGAACGCATCCCCATCGATAACCGAAATTCTTTAATAACAAGACCATGCGGTCTGATTATACCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTTATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA'),
            ('FZTHQMS01DE1KN | cluster size: 2', 'CATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGGTTTGGTGAGCCGTTACCTCACCAACTGCCTAATGGAACGCATCCCCATCGATAACCGAAATTCTTTAATAACAAGACCATGCGGTCTGATTATACCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTTATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA')]

        singletons = [(
            'FZTHQMS01B8T1H',
            'CATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGGTTTGGTGAGCCGTTACCTCACCAACTGCCTAATGGAACGCATCCCCATCGATAACCGAAATTCTTTAATAATTAAACCATGCGGTTTTATTATACCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTTATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCATCACTTA')]

        store_clusters(mapping, self.tiny_test, self.tmpdir)
        actual_centroids = list(
            parse_fasta(open(self.tmpdir + "centroids.fasta")))
        self.assertEqual(actual_centroids, centroids)
        actual_singletons = list(
            parse_fasta(open(self.tmpdir + "singletons.fasta")))
        self.assertEqual(actual_singletons, singletons)

    def test_get_representatives(self):
        """get_representatives should return the representatives as list of Sequence."""

        result = """>1: 5
ABABABA
>3: 1
BABA
>4: 1
ABABAA
>8: 2
BABBA"""
        seqs = self.data.iteritems
        mapping = self.mapping
        test_result = list(get_representatives(mapping, seqs()))
        test_result_as_fasta = "\n".join(
            map(lambda a: a.toFasta(), test_result))

        self.assertEqual(test_result_as_fasta, result)

        # another example
        mapping = {'1': ('a', 'b', 'c'),
                   '2': ('d', 'e', 'f')}
        seqs = [('1', "ACGT"), ('2', "TAGC"), ('a', "TTTTT")]

        observed = list(get_representatives(mapping, seqs))
        expected = [Sequence(name=">1", seq="ACGT"), Sequence(name='2',
                                                              seq="TAGC")]
        self.assertEqual(observed, expected)

    def test_squeeze_seq(self):
        """squeeze should collapse homopolymers to one nuc."""

        seq = "AAAGGGAAACCCGGGA"
        self.assertEqual(squeeze_seq(seq), "AGACGA")
        self.assertEqual(squeeze_seq("AAAATATTTAGGC"), "ATATAGC")
        self.assertEqual(squeeze_seq(""), "")
        self.assertEqual(squeeze_seq("ATGCATGCATGC"), "ATGCATGCATGC")

    def test_wait_for_file(self):
        """wait_for_file should go to sleep if file is not present."""

        # wait_for_file has a debug/test mode, in which it raises an exception
        # instead of going to sleep
        # should not raise anything on valid file
        try:
            wait_for_file("/tmp/denoiser_utils_dummy.tmp", test_mode=True)
        except RuntimeWarning:
            self.fail("wait_for_file fails on valid file")

        # but should raise on file not present
        self.assertRaises(RuntimeWarning, wait_for_file, "/foo/bar/baz",
                          test_mode=True)

    # def test_wait_for_cluster_ids(self):
    #    """wait_for_cluster_ids sleeps until jobs are finished."""
    #
    #   try:
    #       wait_for_cluster_ids([])
    #    except ApplicationNotFoundError:
    #       self.fail("qstat not found. Can't run on cluster.")

        # Can we test a real scenario with submitting a simple sleep script?

    def test_init_flowgram_file(self):
        """init_flowgram_file opens an file and writes header."""
        fh, tmp_filename = init_flowgram_file(n=100, l=400)
        self.assert_(exists(tmp_filename))
        self.tmp_filename = tmp_filename
        fh.close()
        result_file_content = list(open(tmp_filename))

        self.assertEqual(result_file_content, ["100 400\n"])

    def test_append_to_flowgram_file(self):
        """append_to_flowgram_file appends a flowgram to a flowgram file."""

        fh, tmp_filename = init_flowgram_file(n=100, l=400)
        self.assert_(exists(tmp_filename))
        self.tmp_filename = tmp_filename

        flow1 = Flowgram("0 1.2 2.1 3.4 0.02 0.01 1.02 0.08")
        append_to_flowgram_file("test_id", flow1, fh)

        flow2 = Flowgram('0.5 1.0 4.1 0.0 0.0 1.23 0.0 3.1',
                         Name='a', floworder="TACG",
                         header_info={
                             'Bases': 'TACCCCAGGG', 'Clip Qual Right': 7,
                             'Flow Indexes': "1\t2\t3\t3\t3\t3\t6\t8\t8\t8"})
        append_to_flowgram_file("test_id2", flow2, fh, trim=True)
        # close and re-open to read from start, seek might work as well here...
        fh.close()
        fh = open(tmp_filename)
        result_file_content = list(fh)
        self.assertEqual(result_file_content, ["100 400\n",
                                               "test_id 8 0.0 1.2 2.1 3.4 0.02 0.01 1.02 0.08\n",
                                               "test_id2 6 0.5 1.0 4.1 0.0 0.0 1.23\n"])

    def test_cat_sff_files(self):
        "cat_sff_files cats sff_files"""

        expected_bases =  "tcagGCTAACTGTAACCCTCTTGGCACCCACTAAACGCCAATCTTGCTGGAG" +\
            "TGTTTACCAGGCACCCAGCAATGTGAATAGTCActgagcgggctggcaaggc"

       # works with no file
        obs_flows, obs_header = cat_sff_files([])
        self.assertEqual(len(obs_flows), 0)
        self.assertEqual(obs_header, None)

        # works with one file
        obs_flows, obs_header = cat_sff_files([sff_file])
        obs_flows = list(obs_flows)
        self.assertEqual(obs_header['Magic Number'], "0x2E736666")
        self.assertEqual(obs_flows[0].Bases, expected_bases)
        self.assertEqual(len(obs_flows), 2)

        # works with two files
        obs_flows, obs_header = (cat_sff_files([sff_file, sff_file]))
        obs_flows = list(obs_flows)

        self.assertEqual(obs_header['Magic Number'], "0x2E736666")
        self.assertEqual(obs_flows[0].Bases, expected_bases)
        self.assertEqual(obs_flows[2].Bases, expected_bases)
        self.assertEqual(len(obs_flows), 4)

    def test_read_denoiser_mapping(self):
        """read_denoiser_mapping reads correctly"""

        mapping = """1:\t2\t3
4:\t5\t6
7:""".split("\n")
        expected = {'1': ['2', '3'],
                    '4': ['5', '6'],
                    '7': []}
        self.assertEqual(read_denoiser_mapping(mapping),
                         expected)

        # empty mapping gives empty result
        self.assertEqual(read_denoiser_mapping([]), {})

    def test_read_denoiser_mapping_empty_lines(self):
        """read_denoiser_mapping handles empty lines"""

        mapping = """1:\t2\t3
4:\t5\t6
7:
""".split("\n")
        expected = {'1': ['2', '3'],
                    '4': ['5', '6'],
                    '7': []}
        self.assertEqual(read_denoiser_mapping(mapping),
                         expected)

        # empty mapping gives empty result
        self.assertEqual(read_denoiser_mapping([]), {})

    def test_sort_ids(self):
        """sort_ids sorts by abundance"""

        mapping = {"1": ["0", "2", "5", "6"],
                   "3": [],
                   "4": [],
                   "11": [1, 2, 3, 4, 5, 6, 7, 8, 9],
                   "8": ["7"]}

        self.assertEqual(sort_ids(["1", "3", "4", "8", "11"], mapping),
                         ["11", "1", "8", "4", "3"])

    def test_sort_seqs_by_clustersize(self):
        """sort_seqs_by_clustersize works"""

        seqs = {'0': "AAA",
                '1': "AAT",
                '2': "ATT",
                '3': "TTT",
                '4': "TAA",
                '5': "TTA",
                '6': "CCC",
                '7': "GGG",
                '8': "GCG"}

        mapping = {"8": ["7", "6"],
                   "1": ["0", "2", "5"],
                   "4": ["3"]}

        observed = list(sort_seqs_by_clustersize(seqs.iteritems(), mapping))
        expected = [('1', "AAT"), ('8', "GCG"), ('4', 'TAA'), ('7', 'GGG'),
                    ('6', 'CCC'), ('5', 'TTA'), ('3', 'TTT'), ('2', 'ATT'),
                    ('0', 'AAA')]
        self.assertEqual(observed, expected)

    def test_sort_ids(self):
        """sort_ids sorts by abundance"""

        mapping = {"1": ["0", "2", "5", "6"],
                   "3": [],
                   "4": [],
                   "11": [1, 2, 3, 4, 5, 6, 7, 8, 9],
                   "8": ["7"]}

        self.assertEqual(sort_ids(["1", "3", "4", "8", "11"], mapping),
                         ["11", "1", "8", "4", "3"])

    def test_checkpoints(self):
        """storing and loading of checkpoints works"""

        self.tmpdir = mkdtemp(dir="./",
                              suffix="_test_checkpoints/")

        bestscores = dict({1: 0.9,
                           2: 1.1,
                           3: 2.3,
                           4: 99.93232344})

        out_fp = write_checkpoint(
            "Key", 99, self.mapping, [1, 2, 3, 4], bestscores,
            [2, 1, 3, 4],
            self.tmpdir)

        observed = read_checkpoint(out_fp)

        self.assertEqual(observed[0], "Key")
        self.assertEqual(observed[1], 99)
        self.assertEqual(observed[2], self.mapping)
        self.assertEqual(observed[3], [1, 2, 3, 4])
        self.assertEqual(observed[4], bestscores)
        self.assertEqual(observed[5], [2, 1, 3, 4])


class TestFlowgramContainerFile(TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_container(self):
        """FlowgramContainerFile works as expected"""

        fc = FlowgramCollection({'a': '1.0 0.0 0.0 1.0 1.0 1.2 1.2 0.8',
                                 'b': '1.2 1.0 0.0 0.8 1.2 2.4 1.0 0.0'})

        f_container = FlowgramContainerFile(header)

        for f in fc:
            f_container.add(f)

        for f_obs, f_exp in zip(f_container, fc):
            self.assertEqual(str(f_obs), str(f_exp))

        # adding after iter started raises errror
        self.assertRaises(ValueError, f_container.add, f_obs)


class TestFlowgramContainerArray(TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_container(self):
        """FlowgramContainerArray works as expectected"""

        fc = FlowgramCollection({'a': '1.0 0.0 0.0 1.0 1.0 1.2 1.2 0.8',
                                 'b': '1.2 1.0 0.0 0.8 1.2 2.4 1.0 0.0'})

        f_container = FlowgramContainerArray(header)

        for f in fc:
            f_container.add(f)

        for f_obs, f_exp in zip(f_container, fc):
            self.assertEqual(str(f_obs), str(f_exp))


header = {'Version': "0001",
          'Magic Number': '0x2E736666',
          'Index Offset': '7773224',
          'Index Length': '93365',
          '# of Reads': '114',
          'Header Length': '440',
          'Key Length': '4',
          '# of Flows': '400',
          'Flowgram Code': '1',
          'Flow Chars':
          'TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG',
          'Key Sequence': 'TCAG'}


sff_file = """Common Header:
  Magic Number:  0x2E736666
  Version:       0001
  Index Offset:  96099976
  Index Length:  1158685
  # of Reads:    57902
  Header Length: 440
  Key Length:    4
  # of Flows:    400
  Flowgram Code: 1
  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
  Key Sequence:  TCAG

>FIQU8OX05GCVRO
  Run Prefix:   R_2008_10_15_16_11_02_
  Region #:     5
  XY Location:  2489_3906

  Run Name:       R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford
  Analysis Name:  /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis
  Full Path:      /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis

  Read Header Len:  32
  Name Length:      14
  # of Bases:       104
  Clip Qual Left:   5
  Clip Qual Right:  85
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.06	0.08	1.04	0.08	0.05	0.94	0.10	2.01	0.10	0.07	0.96	0.09	1.04	1.96	1.07	0.10	1.01	0.13	0.08	1.01	1.06	1.83	2.89	0.18	0.96	0.13	0.99	0.11	1.94	0.12	0.13	1.92	0.21	0.07	0.94	0.17	0.03	0.97	2.76	0.15	0.05	1.02	1.14	0.10	0.98	2.54	1.13	0.96	0.15	0.21	1.90	0.16	0.07	1.78	0.22	0.07	0.93	0.22	0.97	0.08	2.02	0.15	0.19	1.02	0.19	0.09	1.02	0.17	0.99	0.09	0.18	1.84	0.16	0.91	0.10	1.10	1.00	0.20	0.09	1.11	3.01	1.07	1.98	0.14	0.22	1.09	0.17	1.99	0.15	0.20	0.92	0.17	0.07	1.01	2.96	0.15	0.07	1.06	0.20	1.00	0.10	0.12	1.00	0.15	0.08	1.90	0.19	0.10	0.99	0.18	0.09	0.99	1.08	0.15	0.07	1.06	0.14	1.84	0.13	0.11	0.95	1.05	0.13	1.04	1.10	0.18	0.94	0.14	0.10	0.97	1.08	0.12	1.08	0.18	0.08	1.00	0.13	0.98	0.15	0.87	0.13	0.19	1.01	3.06	0.17	0.11	1.04	0.09	1.03	0.10	0.11	2.02	0.16	0.11	1.04	0.04	0.09	1.87	0.13	2.09	0.13	0.10	0.97	0.17	0.08	0.08	0.04	0.12	0.05	0.08	0.07	0.08	0.05	0.07	0.06	0.07	0.03	0.05	0.04	0.09	0.04	0.07	0.04	0.07	0.06	0.03	0.06	0.06	0.06	0.06	0.07	0.09	0.04	0.05	0.08	0.05	0.04	0.09	0.06	0.03	0.02	0.08	0.04	0.06	0.05	0.08	0.03	0.08	0.05	0.05	0.05	0.10	0.05	0.05	0.07	0.06	0.04	0.06	0.05	0.03	0.04	0.05	0.06	0.04	0.04	0.07	0.04	0.04	0.05	0.05	0.04	0.07	0.06	0.05	0.03	0.08	0.05	0.06	0.04	0.06	0.05	0.04	0.04	0.04	0.05	0.06	0.04	0.05	0.04	0.05	0.05	0.06	0.05	0.06	0.04	0.06	0.07	0.06	0.05	0.05	0.05	0.06	0.06	0.04	0.05	0.06	0.03	0.06	0.04	0.06	0.05	0.03	0.06	0.06	0.05	0.06	0.04	0.03	0.06	0.06	0.06	0.03	0.04	0.05	0.05	0.07	0.04	0.05	0.06	0.07	0.07	0.05	0.07	0.06	0.05	0.06	0.05	0.07	0.06	0.05	0.06	0.07	0.05	0.06	0.04	0.06	0.05	0.05	0.06	0.04	0.06	0.04	0.03	0.06	0.05	0.05	0.04	0.05	0.05	0.04	0.04	0.05	0.06	0.06	0.04	0.04	0.05	0.06	0.04	0.04	0.04	0.05	0.05	0.04	0.05	0.05	0.03	0.06	0.06	0.06	0.04	0.07	0.05	0.05	0.04	0.06	0.06	0.05	0.05	0.07	0.04	0.06	0.06	0.06	0.04	0.06	0.03	0.06	0.04	0.06	0.04	0.09	0.05	0.05	0.05	0.07	0.06	0.05	0.05	0.06	0.05	0.05	0.05	0.04	0.04	0.06	0.05	0.05	0.05	0.05	0.04	0.05	0.05	0.06	0.04	0.05	0.05	0.05	0.05	0.05	0.04	0.06	0.04	0.05	0.05	0.04	0.05	0.05	0.05	0.04
Flow Indexes:	1	3	6	8	8	11	13	14	14	15	17	20	21	22	22	23	23	23	25	27	29	29	32	32	35	38	39	39	39	42	43	45	46	46	46	47	48	51	51	54	54	57	59	61	61	64	67	69	72	72	74	76	77	80	81	81	81	82	83	83	86	88	88	91	94	95	95	95	98	100	103	106	106	109	112	113	116	118	118	121	122	124	125	127	130	131	133	136	138	140	143	144	144	144	147	149	152	152	155	158	158	160	160	163
Bases:	tcagGCTAACTGTAACCCTCTTGGCACCCACTAAACGCCAATCTTGCTGGAGTGTTTACCAGGCACCCAGCAATGTGAATAGTCActgagcgggctggcaaggc
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	37	37	37	37	37	39	39	39	39	24	24	24	37	34	28	24	24	24	28	34	39	39	39	39	39	39	39	39	39	39	39	39	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37

>FIQU8OX05F8ILF
  Run Prefix:   R_2008_10_15_16_11_02_
  Region #:     5
  XY Location:  2440_0913

  Run Name:       R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford
  Analysis Name:  /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis
  Full Path:      /data/2008_10_15/R_2008_10_15_16_11_02_FLX04070166_adminrig_1548jinnescurtisstanford/D_2008_10_15_15_12_26_FLX04070166_1548jinnescurtisstanford_FullAnalysis

  Read Header Len:  32
  Name Length:      14
  # of Bases:       206
  Clip Qual Left:   5
  Clip Qual Right:  187
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.00	1.01	0.00	0.00	1.00	0.00	1.00	0.00	1.05	0.00	0.91	0.10	1.07	0.95	1.01	0.00	0.06	0.93	0.02	0.03	1.06	1.18	0.09	1.00	0.05	0.90	0.11	0.07	1.99	0.11	0.02	1.96	1.04	0.13	0.01	2.83	0.10	1.97	0.06	0.11	1.04	0.13	0.03	0.98	1.15	0.07	1.00	0.07	0.08	0.98	0.11	1.92	0.05	0.04	2.96	1.02	1.02	0.04	0.93	1.00	0.13	0.04	1.00	1.03	0.08	0.97	0.13	0.11	1.88	0.09	0.05	1.02	1.89	0.07	0.11	0.98	0.05	0.07	1.01	0.08	0.05	1.01	0.13	1.00	0.07	0.10	1.04	0.10	0.04	0.98	0.12	1.03	0.96	0.11	0.07	1.00	0.09	0.03	1.03	0.11	1.95	1.06	0.13	0.05	1.00	0.13	0.11	1.00	0.09	0.03	2.89	0.08	0.95	0.09	1.03	1.02	1.05	1.07	0.08	0.12	2.81	0.08	0.08	1.00	1.07	0.07	0.05	1.86	0.12	0.98	0.06	2.00	0.11	1.02	0.11	0.08	1.88	0.13	1.03	0.13	0.98	0.15	0.11	1.03	1.03	1.04	0.18	0.98	0.13	0.15	1.04	0.11	1.01	0.13	0.06	1.01	0.06	1.02	0.08	0.99	0.14	0.99	0.09	0.05	1.09	0.04	0.07	2.96	0.09	2.03	0.13	2.96	1.13	0.08	1.03	0.07	0.99	0.11	0.05	1.05	1.04	0.09	0.07	1.00	1.03	0.09	0.06	1.06	1.04	2.94	0.18	0.06	0.93	0.10	1.10	0.11	2.02	0.17	1.00	1.03	0.06	0.11	0.96	0.04	3.00	0.11	0.07	1.99	0.10	2.03	0.12	0.97	0.16	0.01	2.09	0.14	1.04	0.16	0.06	1.03	0.14	1.12	0.12	0.05	0.96	1.01	0.10	0.14	0.94	0.03	0.12	1.10	0.92	0.09	1.10	1.04	1.02	0.12	0.97	2.00	0.15	1.08	0.04	1.03	1.04	0.03	0.09	5.16	1.02	0.09	0.13	2.66	0.09	0.05	1.06	0.07	0.89	0.05	0.12	1.10	0.16	0.06	1.01	0.13	1.00	0.14	0.98	0.09	2.92	1.28	0.03	2.95	0.98	0.16	0.08	0.95	0.96	1.09	0.08	1.07	1.01	0.16	0.06	4.52	0.12	1.03	0.07	0.09	1.03	0.14	0.03	1.01	1.99	1.05	0.14	1.03	0.13	0.03	1.10	0.10	0.96	0.11	0.99	0.12	0.05	0.94	2.83	0.14	0.12	0.96	0.00	1.00	0.11	0.14	1.98	0.08	0.11	1.04	0.01	0.11	2.03	0.15	2.05	0.10	0.03	0.93	0.01	0.08	0.12	0.00	0.16	0.05	0.07	0.08	0.11	0.07	0.05	0.04	0.10	0.05	0.05	0.03	0.07	0.03	0.04	0.04	0.06	0.03	0.05	0.04	0.09	0.03	0.08	0.03	0.07	0.02	0.05	0.02	0.06	0.01	0.05	0.04	0.06	0.02	0.04	0.04	0.04	0.03	0.03	0.06	0.06	0.03	0.02	0.02	0.08	0.03	0.01	0.01	0.06	0.03	0.01	0.03	0.04	0.02	0.00	0.02	0.05	0.00	0.02	0.02	0.03	0.00	0.02	0.02	0.04	0.01	0.00	0.01	0.05
Flow Indexes:	1	3	6	8	10	12	14	15	16	19	22	23	25	27	30	30	33	33	34	37	37	37	39	39	42	45	46	48	51	53	53	56	56	56	57	58	60	61	64	65	67	70	70	73	74	74	77	80	83	85	88	91	93	94	97	100	102	102	103	106	109	112	112	112	114	116	117	118	119	122	122	122	125	126	129	129	131	133	133	135	138	138	140	142	145	146	147	149	152	154	157	159	161	163	166	169	169	169	171	171	173	173	173	174	176	178	181	182	185	186	189	190	191	191	191	194	196	198	198	200	201	204	206	206	206	209	209	211	211	213	216	216	218	221	223	226	227	230	233	234	236	237	238	240	241	241	243	245	246	249	249	249	249	249	250	253	253	253	256	258	261	264	266	268	270	270	270	271	273	273	273	274	277	278	279	281	282	285	285	285	285	285	287	290	293	294	294	295	297	300	302	304	307	308	308	308	311	313	316	316	319	322	322	324	324	327
Bases:	tcagAGACGCACTCAATTATTTCCATAGCTTGGGTAGTGTCAATAATGCTGCTATGAACATGGGAGTACAAATATTCTTCAAGATACTGATCTCATTTCCTTTAGATATATACCCAGAAGTGAAATTCCTGGATCACATAGTAGTTCTATTTTTATTTGATGAGAAACTTTATACTATTTTTCATAActgagcgggctggcaaggc
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	34	34	34	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	36	36	38	25	25	25	38	37	37	37	37	37	37	33	33	34	37	37	37	37	37	37	37	38	34	20	20	26	26	20	34	38	37	37	37	37	37	37	37	37	37	38	38	38	37	37	37	37	37	37	37	37	37	37

""".split('\n')

if __name__ == "__main__":
    main()
