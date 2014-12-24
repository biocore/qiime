#!/usr/bin/env python
# file test_exclude_seqs_by_blast.py

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Zaneveld", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.9.0-rc1"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"

"""
Test code for exclude_seqs_by_blast.py.

NOTE: requires BLAST to be properly installed with
environment variable set for tests to pass
"""
from os import remove, system, close
from random import choice
from tempfile import mkstemp

from numpy import array, arange, log, log10
from unittest import TestCase, main
from numpy.testing import assert_almost_equal
from bfillings.blast import BlastResult
from qiime.exclude_seqs_by_blast import blast_genome,\
    find_homologs,\
    sequences_to_file,\
    no_filter,\
    make_percent_align_filter,\
    query_ids_from_blast_result,\
    ids_from_fasta_lines,\
    id_from_fasta_label_line,\
    seqs_from_file,\
    ids_to_seq_file
from qiime.util import remove_files


class ExcludeHumanTests(TestCase):

    def setUp(self):

        self.blast_lines = BLAST_LINES
        self.blast_result = BlastResult(self.blast_lines)

        fd, self.subjectdb_fp = mkstemp(prefix='ExcludeByBlastTests_',
                                       suffix='.fasta')
        close(fd)
        fd, self.query_fp = mkstemp(prefix='ExcludeByBlastTests_',
                                   suffix='.fasta')
        close(fd)
        fd, self.query2_fp = mkstemp(prefix='ExcludeByBlastTests_',
                                    suffix='.fasta')
        close(fd)

        open(self.subjectdb_fp, "w").writelines(TEST_BLAST_DB_LINES)
        open(self.query_fp, "w").writelines(TEST_BLAST_DB_LINES)
        open(self.query2_fp, "w").writelines(TEST_BLAST_DB2_LINES)

        self._paths_to_clean_up = [self.subjectdb_fp, self.query_fp,
                                   self.query2_fp]

    def tearDown(self):
        remove_files(self._paths_to_clean_up)

    def test_blast_genome(self):
        """blast_genome should return raw BLAST output."""

        formatdb_cmd = 'formatdb -p F -o T -i %s' % self.subjectdb_fp
        system(formatdb_cmd)
        self._paths_to_clean_up.append("formatdb.log")
        for suffix in ["nhr", "nin", "nsd", "nsi", "nsq"]:
            self._paths_to_clean_up.append(".".join(
                [self.subjectdb_fp, suffix]))

        raw_output = blast_genome(TEST_BLAST_DB_LINES, self.subjectdb_fp,
                                  e_value=1e-4, max_hits=100, word_size=28,
                                  working_dir="./", blast_mat_root=None)

        i = 0
        for line in raw_output:

            if line.startswith("#"):
                i += 1
                continue  # comments depend on tmpfilename, BLAST version
            self.assertEqual(raw_output[i], EXP_BLAST_OUTPUT[i])
            i += 1

    def test_find_homologs(self):
        """find_homologs should return raw data, filtered and removed ids."""

        formatdb_cmd = 'formatdb -p F -o T -i %s' % self.subjectdb_fp
        system(formatdb_cmd)
        self._paths_to_clean_up.append("formatdb.log")
        for suffix in ["nhr", "nin", "nsd", "nsi", "nsq"]:
            self._paths_to_clean_up.append(".".join(
                [self.subjectdb_fp, suffix]))

        blast_output, hit_ids, removed_hit_ids =\
            find_homologs(self.query_fp, self.subjectdb_fp, e_value=1e-4,
                          max_hits=100, working_dir="./", blast_mat_root=None,
                          wordsize=28, percent_aligned=0.98, DEBUG=False)

        self.assertEqual(hit_ids, set(["bth:BT_0001", "hsa:8355"]))
        self.assertEqual(removed_hit_ids, set())

        i = 0
        for line in blast_output:

            if line.startswith("#"):
                i += 1
                continue  # depends on tmpfilename, skip testing

            self.assertEqual(blast_output[i], EXP_BLAST_OUTPUT[i])
            i += 1

        # Ensure low % alignment seqs are removed
        blast_output, hit_ids, removed_hit_ids =\
            find_homologs(self.query2_fp, self.subjectdb_fp,
                          e_value=1e-4, max_hits=100, working_dir="./",
                          blast_mat_root=None, wordsize=28, percent_aligned=1.00,
                          DEBUG=False)

        self.assertEqual(hit_ids, set(["bth:BT_0001"]))
        self.assertEqual(removed_hit_ids, set(["hsa:8355_tweaked"]))

        # Ensure high % alignment seqs are not removed
        blast_output, hit_ids, removed_hit_ids =\
            find_homologs(self.query2_fp, self.subjectdb_fp,
                          e_value=1e-4, max_hits=100, working_dir="./",
                          blast_mat_root=None, wordsize=28, percent_aligned=0.75,
                          DEBUG=False)

        self.assertEqual(hit_ids, set(["bth:BT_0001", "hsa:8355_tweaked"]))
        self.assertEqual(removed_hit_ids, set())

    def test_sequences_to_file(self):
        """sequences_to_file should write a standard format FASTA file."""

        fd, self.seq_test_fp = mkstemp(prefix='ExcludeByBlastTests_',
                                      suffix='.fasta')
        close(fd)
        self._paths_to_clean_up.append(self.seq_test_fp)

        ids = ["bth:BT_0001", "hsa:8355"]
        seqs = seqs_from_file(ids, open(self.query_fp).readlines())
        sequences_to_file(seqs, self.seq_test_fp)

        self.assertEqual(open(self.seq_test_fp).readlines(),
                         open(self.query_fp).readlines())

    def test_no_filter(self):
        """no_filter should always return True."""

        d1 = {"% IDENTITY": "97.6"}
        d2 = {"% IDENTITY": "0.0"}
        d3 = {"% IDENTITY": "100.0"}

        self.assertTrue(no_filter(d1))
        self.assertTrue(no_filter(d2))
        self.assertTrue(no_filter(d3))

    def test_make_percent_align_filter(self):
        """make_percent_align_filter should return a percent align filter fn"""

        d1 = {"% IDENTITY": "97.6"}
        d2 = {"% IDENTITY": "0.0"}
        d3 = {"% IDENTITY": "100.0"}

        af1 = make_percent_align_filter(0.50)
        af2 = make_percent_align_filter(0.00)
        af3 = make_percent_align_filter(1.0)

        # Test filter 1
        self.assertTrue(af1(d1))
        self.assertFalse(af1(d2))
        self.assertTrue(af1(d3))

        # Test filter 2
        self.assertTrue(af2(d1))
        self.assertTrue(af2(d2))
        self.assertTrue(af2(d3))

        # Test filter 3
        self.assertFalse(af3(d1))
        self.assertFalse(af3(d2))
        self.assertTrue(af3(d3))

    def test_query_ids_from_blast_result(self):
        "query_ids_from_blast_result should return query_ids matching filter"
        align_filter = make_percent_align_filter(2.0)  # none should pass

        ok_ids, removed_ids = query_ids_from_blast_result(
            self.blast_result, align_filter, DEBUG=True)
        self.assertEqual(ok_ids, set())

    def test_ids_from_fasta_lines(self):
        """ ids_from_fasta_lines should return ids"""

        fasta_lines = \
            [">hsa:8355  HIST1H3G; histone cluster 1, H3g ; K11253 histone H3",
             "atggcccgcaccaagcagactgcacgcaagtccaccggtggcaaagcgccgcgcaagcagctgg",
             "ccactaaggcggctcggaaaagcgcgccggccaccggcggcgtgaagaaacctcatcgctaccg",
             "tcccggcaccgtggctctgcgcgagattcgccgctatcagaagtcgactgagctgctgatccgc",
             "aagttgcctttccaacgcctggtgcgagaaatcgctcaggacttcaagacagatctgcgctttc",
             "agagttccgcggtgatggccctgcaggaggcctgcgaggcctacttggtggggctctttgagga",
             "taccaacctgtgtgccatccatgctaagcgagtgactatcatgcccaaggacattcagctcgct",
             "cgccgcattcgtggggagagagcgtag",
             ">hsa:9081  PRY; PTPN13-like, Y-linked",
             "atgggagccactgggcttggctttctactttcctggagacaagacaatttgaatggcact"]
        exp_ids = ["hsa:8355", "hsa:9081"]
        obs_ids = ids_from_fasta_lines(fasta_lines)

        self.assertEqual(obs_ids, exp_ids)

    def test_id_from_fasta_label_line(self):
        """id_from_fasta_label_line should extract id"""
        label_line = \
            ">hsa:8355  HIST1H3G; histone cluster 1, H3g ; K11253 histone H3"
        self.assertEqual(id_from_fasta_label_line(label_line), "hsa:8355")

    def test_seqs_from_file(self):
        """seqs_from_file should extract labels,seqs for specified ids"""
        ids = "bth:BT_0001"
        seqs = seqs_from_file(ids, open(self.query_fp).readlines())
        all_results = []
        for label, seq in seqs:
            all_results.append((label, seq))

        self.assertEqual(1, len(all_results))  # should return only 1 entry

        # fn should return version lacking ">" and newlines
        self.assertEqual(all_results[0],
                         (TEST_BLAST_DB_LINES[2].strip(">").strip(),
                          TEST_BLAST_DB_LINES[3].strip()))

    def test_ids_to_seq_file(self):
        """ids_to_seq_file should lookup and write out seqs for given ids"""
        fd, self.id_test_fp = mkstemp(prefix='ExcludeByBlastTests_',
                                     suffix='.fasta')
        close(fd)

        self._paths_to_clean_up.append(self.id_test_fp)

        ids = ["bth:BT_0001"]
        ids_to_seq_file(ids, self.query_fp, self.id_test_fp)

        # this is the bth entry
        exp_lines = open(self.query_fp).readlines()[2:]
        self.assertEqual(open(self.id_test_fp).readlines(), exp_lines)


# Predefined Test strings
BLAST_LINES = [
    '# BLASTN 2.2.16 [Mar-25-2007]\n',
    '# Query: hsa:8355  HIST1H3G; histone cluster 1, H3g ; K11253 histone H3\n',
    '# Database: /home/zaneveld/quicksand/data/human_genome/h.sapiens.nuc\n',
    '# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n',
    'hsa:8355\thsa:8355\t95.00\t411\t0\t0\t1\t411\t1\t411\t0.0\t 815\n',
    'hsa:8355\thsa:8351\t88.29\t410\t48\t0\t1\t410\t1\t410\t8e-121\t 432\n',
    'hsa:8355\thsa:8353\t87.15\t397\t51\t0\t13\t409\t13\t409\t7e-106\t 383\n',
    'hsa:8355\thsa:8968\t86.63\t404\t54\t0\t1\t404\t1\t404\t6e-103\t 373\n',
    'hsa:8355\thsa:8354\t86.84\t380\t50\t0\t25\t404\t25\t404\t4e-98\t 357\n',
    '# BLASTN 2.2.16 [Mar-25-2007]\n',
    '# Query: hsa:9081  PRY; PTPN13-like, Y-linked\n',
    '# Database: /home/zaneveld/quicksand/data/human_genome/h.sapiens.nuc\n',
    '# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n',
    'hsa:9081\thsa:442862\t100.00\t444\t0\t0\t1\t444\t1\t444\t0.0\t 880\n',
    'hsa:9081\thsa:9081\t100.00\t444\t0\t0\t1\t444\t1\t444\t0.0\t 880\n',
    '# BLASTN 2.2.16 [Mar-25-2007]\n',
    '# Query: hsa:23434  C3orf27; chromosome 3 open reading frame 27\n',
    '# Database: /home/zaneveld/quicksand/data/human_genome/h.sapiens.nuc\n',
    '# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n',
    'hsa:23434\thsa:23434\t100.00\t236\t0\t0\t1\t236\t1\t236\t2e-131\t 468\n',
    'hsa:23434\thsa:23434\t100.00\t196\t0\t0\t255\t450\t255\t450\t1e-107\t 389\n']

TEST_BLAST_DB2_LINES = [
    ">hsa:8355_tweaked HIST1H3G; <fake data for testing>\n",
    """atggcccgcaccaagcagactgcacgcaagtccaccggtggcaaagcgccgcgcaagcagctggccactaaggcggctttggaaaagcgcgccggccaccggcggcgtgaagaaacctcatcgctaccgtcccggcaccgtggctctgcgcgagattcgccgctatcagaagtcgactgagctgctgatccgcaagttgcctttccaacgcctggtgcgagaaatcgctcaggacttcaagacagatctgcgctttcagacttccgcggtgatggccctgcaggaggcctgcgaggcctacttggtggggctctttgaggataccaacctgtgtgccatccatgctaagcgagtgactatcatgcccaaggacattcagctcgctcgccgcattcgtggggagagagcgtag\n""",
    ">bth:BT_0001  hypothetical protein\n",
    """ttggtatctaccagtacgcacgacgatgcttttgacttcgactttggttacactggtaagcttcagttcttggtagccactgtagatgcaaatagtacctattacactaaagacccgaatggtattgaatgtgataacgacggaagcagttcatctttaactccgttcactcacccgacaatcagtaacttaacaatcgttggaaccgttaatggtaaggttgcacaatctgcaatgggtgatggtaaatccatgaaatcttgtgccaacttccgtagaaactgccaatttactttggtgaacagtattctttacggatatcctaccggtatcttgtgtgaaaccactaacagctatgttttcaaaaacaatgttgtaaatggtgttagtactacattttcaggtatcacagctgacgcgactaatactgctgctgcaagtgctgaggctattgggctgacttctccgtggggtggatatacaggtttgatgcctaatgcatctccagccaatgcaggtgcagattttagtgaattggatagttggtttacgactacttcttacagaggtgctgttggtggacgttcaaactggttaactcaagcgtgggtaaaataa\n"""]


TEST_BLAST_DB_LINES = [
    ">hsa:8355  HIST1H3G; histone cluster 1, H3g ; K11253 histone H3\n",
    """atggcccgcaccaagcagactgcacgcaagtccaccggtggcaaagcgccgcgcaagcagctggccactaaggcggctcggaaaagcgcgccggccaccggcggcgtgaagaaacctcatcgctaccgtcccggcaccgtggctctgcgcgagattcgccgctatcagaagtcgactgagctgctgatccgcaagttgcctttccaacgcctggtgcgagaaatcgctcaggacttcaagacagatctgcgctttcagagttccgcggtgatggccctgcaggaggcctgcgaggcctacttggtggggctctttgaggataccaacctgtgtgccatccatgctaagcgagtgactatcatgcccaaggacattcagctcgctcgccgcattcgtggggagagagcgtag\n""",
    ">bth:BT_0001  hypothetical protein\n",
    """ttggtatctaccagtacgcacgacgatgcttttgacttcgactttggttacactggtaagcttcagttcttggtagccactgtagatgcaaatagtacctattacactaaagacccgaatggtattgaatgtgataacgacggaagcagttcatctttaactccgttcactcacccgacaatcagtaacttaacaatcgttggaaccgttaatggtaaggttgcacaatctgcaatgggtgatggtaaatccatgaaatcttgtgccaacttccgtagaaactgccaatttactttggtgaacagtattctttacggatatcctaccggtatcttgtgtgaaaccactaacagctatgttttcaaaaacaatgttgtaaatggtgttagtactacattttcaggtatcacagctgacgcgactaatactgctgctgcaagtgctgaggctattgggctgacttctccgtggggtggatatacaggtttgatgcctaatgcatctccagccaatgcaggtgcagattttagtgaattggatagttggtttacgactacttcttacagaggtgctgttggtggacgttcaaactggttaactcaagcgtgggtaaaataa\n"""]

EXP_BLAST_OUTPUT =\
    ['# BLASTN 2.2.16 [Mar-25-2007]\n',
     '# Query: hsa:8355  HIST1H3G; histone cluster 1, H3g ; K11253 histone H3\n',
     '# Database: /tmp/ExcludeByBlastTests_FQifdUaBoAngS1XGIvnP.fasta\n',
     '# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n',
     'hsa:8355\thsa:8355\t100.00\t411\t0\t0\t1\t411\t1\t411\t0.0\t 815\n',
     '# BLASTN 2.2.16 [Mar-25-2007]\n',
     '# Query: bth:BT_0001  hypothetical protein\n',
     '# Database: /tmp/ExcludeByBlastTests_FQifdUaBoAngS1XGIvnP.fasta\n',
     '# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n',
     'bth:BT_0001\tbth:BT_0001\t100.00\t618\t0\t0\t1\t618\t1\t618\t0.0\t1225\n']

if __name__ == "__main__":
    main()
