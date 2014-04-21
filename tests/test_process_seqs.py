#!/usr/bin/env python

from unittest import TestCase, main

import numpy as np

from future.builtins import zip

from skbio.core.iterator import FastqIterator
from skbio.parse.sequences import parse_fastq
from qiime.process_seqs import IterAdapter, SequenceWorkflow, count_mismatches
from qiime.util import MetadataMap


class SupportTests(TestCase):
    def setUp(self):
        pass

    def test_count_mismatches(self):
        """Count mismatches in sequences"""
        s1 = "AATTGGCC"
        s2 = "AATTCCCC"

        self.assertEqual(count_mismatches(s1, s1), 0)
        self.assertEqual(count_mismatches(s1, s2), 2)
        self.assertEqual(count_mismatches(s2, s1), 2)
        self.assertEqual(count_mismatches(s2, s2), 0)

class IterAdapterTests(TestCase):
    def test_iter(self):
        seq_raw = fastq1.splitlines()
        bc_raw = barcode_fastq1.splitlines()

        seq = FastqIterator([seq_raw])
        barcode = FastqIterator([bc_raw])
        it = IterAdapter(seq=seq, barcode=barcode)

        for rec, s, b in zip(it, parse_fastq(seq_raw), parse_fastq(bc_raw)):
            self.assertEqual(rec['SequenceID'], s[0])
            self.assertEqual(rec['Sequence'], s[1])
            np.testing.assert_equal(rec['Qual'], s[2])
            self.assertEqual(rec['BarcodeID'], b[0])
            self.assertEqual(rec['Barcode'], b[1])
            np.testing.assert_equal(rec['BarcodeQual'], b[2])


class ProcessSeqsWorkflowTests(TestCase):
    """Basing structure off of test_split_libraries_fastq.py"""
    def setUp(self):
        self.fastq1 = fastq1.split('\n')
        self.barcode_fastq1 = barcode_fastq1.split('\n')
        self.fastq2 = fastq2.split('\n')
        self.barcode_fastq2 = barcode_fastq2.split('\n')
        self.fastq1_expected_no_qual_unassigned = fastq1_expected_no_qual_unassigned
        self.fastq1_expected_default = fastq1_expected_default
        self.fastq2_expected_default = fastq2_expected_default
        self.fastq1_expected_single_barcode = fastq1_expected_single_barcode
        self.mapping = mapping

    def _make_workflow_obj(self, options):
        """Helper method for creating workflows"""
        return SequenceWorkflow(Options=options, Mapping=self.mapping)

    def test_workflow_construction(self):
        """Make sure we can construct using our helper method"""
        x = self._make_workflow_obj({'foo':'bar'})

    def test_wf_init(self):
        """Check the initialization method"""
        wf_obj = self._make_workflow_obj({'foo':'bar'})
        wf_obj.FinalState['Sequence'] = 'w00t'
        wf_obj.wf_init({'Sequence':'foo'})
        self.assertEqual(set(wf_obj.FinalState.values()), set([None, 'foo']))

    def test_quality_max_bad_run_length(self):
        """Verify max bad run length quality trimming"""
        wf_obj = self._make_workflow_obj({'phred_quality_threshold':5,
                                          'max_bad_run_length':3})
        item1 = {'Sequence':'AATTGGCC',
                 'Qual':array([6, 6, 6, 6, 6, 6, 6, 6])}
        exp1 = item1.copy()

        item2 = {'Sequence':'AATTGGCC',
                 'Qual':array([6, 6, 6, 1, 1, 6, 6, 6])}
        exp2 = item2.copy()

        item3 = {'Sequence':'AATTGGCC',
                 'Qual':array([6, 6, 1, 1, 1, 1, 6, 6])}
        exp3 = {'Sequence':'AA', 'Qual':array([6, 6])}

        wf_obj._quality_max_bad_run_length(item1)
        wf_obj._quality_max_bad_run_length(item2)
        wf_obj._quality_max_bad_run_length(item3)

        self.assertEqual(item1, exp1)
        self.assertEqual(item2, exp2)
        self.assertEqual(item3, exp3)

    def test_quality_min_per_read_length_fraction(self):
        """Verify minimum quality per read length"""
        wf_obj = self._make_workflow_obj({'phred_quality_threshold':5,
                                          'min_per_read_length_fraction':0.6})
        item1 = {'Sequence':'AATTGGCC',
                 'Qual':array([6, 6, 6, 6, 6, 6, 6, 6])}
        exp1 = item1.copy()

        item2 = {'Sequence':'AATTGGCC',
                 'Qual':array([6, 1, 6, 1, 1, 6, 6, 6])}
        exp2 = item2.copy()

        item3 = {'Sequence':'AATTGGCC',
                 'Qual':array([6, 6, 1, 1, 1, 1, 6, 6])}
        exp3 = {'Sequence':'AATTGGCC', 'Qual':array([6, 6, 1, 1, 1, 1, 6, 6])}

        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._quality_min_per_read_length_fraction(item1)
        self.assertFalse(wf_obj.Failed)

        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._quality_min_per_read_length_fraction(item2)
        self.assertFalse(wf_obj.Failed)

        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._quality_min_per_read_length_fraction(item3)
        self.assertTrue(wf_obj.Failed)

        self.assertEqual(item1, exp1)
        self.assertEqual(item2, exp2)
        self.assertEqual(item3, exp3)

    def test_demultiplex_golay12(self):
        # this is a wrapper, tested in test_deultiplex_encoded_barcode
        pass

    def test_demultiplex_hamming8(self):
        # this is a wrapper, tested in test_deultiplex_encoded_barcode
        pass

    def test_demultiplex_encoded_barcode(self):
        """Verify decoding barcodes"""
        wf_obj = self._make_workflow_obj({})

        needs_a_fix = {'Barcode':'GGAGACAAGGGT', 'Sequence':'AATTGGCC'}
        exact = {'Barcode':'GGAGACAAGGGA', 'Sequence':'AATTGGCC'}
        from_sequence = {'Barcode':None, 'Sequence':'GGAGACAAGGGAAATTAATT'}
        unknown_barcode = {'Barcode':'ACACCTGGTGAT', 'Sequence':'AATTGGCC'}

        wf_obj.wf_init(needs_a_fix)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._demultiplex_encoded_barcode(needs_a_fix)
        self.assertEqual(wf_obj.FinalState['Original barcode'], 'GGAGACAAGGGT')
        self.assertEqual(wf_obj.FinalState['Corrected barcode errors'], 1)
        self.assertEqual(wf_obj.FinalState['Corrected barcode'], 'GGAGACAAGGGA')
        self.assertEqual(wf_obj.FinalState['Sample'], 's5')
        self.assertFalse(wf_obj.Failed)

        wf_obj.wf_init(exact)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._demultiplex_encoded_barcode(exact)
        self.assertEqual(wf_obj.FinalState['Original barcode'], 'GGAGACAAGGGA')
        self.assertEqual(wf_obj.FinalState['Corrected barcode errors'], 0)
        self.assertEqual(wf_obj.FinalState['Corrected barcode'], None)
        self.assertEqual(wf_obj.FinalState['Sample'], 's5')
        self.assertFalse(wf_obj.Failed)

        wf_obj.wf_init(from_sequence)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._demultiplex_encoded_barcode(from_sequence)
        self.assertEqual(wf_obj.FinalState['Original barcode'], 'GGAGACAAGGGA')
        self.assertEqual(wf_obj.FinalState['Corrected barcode errors'], 0)
        self.assertEqual(wf_obj.FinalState['Corrected barcode'], None)
        self.assertEqual(wf_obj.FinalState['Sample'], 's5')
        self.assertFalse(wf_obj.Failed)

        wf_obj.wf_init(unknown_barcode)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._demultiplex_encoded_barcode(unknown_barcode)
        self.assertEqual(wf_obj.FinalState['Original barcode'], 'ACACCTGGTGAT')
        self.assertEqual(wf_obj.FinalState['Corrected barcode errors'], 0)
        self.assertEqual(wf_obj.FinalState['Corrected barcode'], 'ACACCTGGTGAT')
        self.assertEqual(wf_obj.FinalState['Sample'], None)
        self.assertTrue(wf_obj.Failed)

    def test_demultiplex_max_barcode_error(self):
        """Verify failing max_barcode_error checking"""
        wf_obj = self._make_workflow_obj({'max_barcode_error':0})

        needs_a_fix = {'Barcode':'GGAGACAAGGGT', 'Sequence':'AATTGGCC'}
        exact = {'Barcode':'GGAGACAAGGGA', 'Sequence':'AATTGGCC'}

        wf_obj.wf_init(exact)
        wf_obj._demultiplex_encoded_barcode(exact)
        wf_obj._demultiplex_max_barcode_error(exact)
        self.assertFalse(wf_obj.Failed)

        wf_obj.wf_init(needs_a_fix)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._demultiplex_encoded_barcode(needs_a_fix)
        self.assertFalse(wf_obj.Failed)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._demultiplex_max_barcode_error(needs_a_fix)
        self.assertTrue(wf_obj.Failed)

    def test_primer_instrument_454(self):
        # individual tests for each method call by this function
        pass

    def test_primer_check_forward(self):
        """Pull the forward primer as expected"""
        wf_obj = self._make_workflow_obj({'max_primer_mismatch':2,
                                          'retain_primer':False})
        item1 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'AATTGGCC',
                 'Qual':array([1,2,3,4,5,6,7,8])}
        item2 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'AATTGCCC',
                 'Qual':array([1,2,3,4,5,6,7,8])}
        item3 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'GGTTGCCC',
                 'Qual':array([1,2,3,4,5,6,7,8])}
        exp_item1 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'CC',
                 'Qual':array([7,8])}
        exp_item2 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'CC',
                 'Qual':array([7,8])}
        exp_item3 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'GGTTGCCC',
                 'Qual':array([1,2,3,4,5,6,7,8])}

        # item is modified in place in these operations as retain_primer is False
        wf_obj.wf_init(item1)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj.FinalState['Final barcode'] = 'AAAAAAAAAAAA'
        wf_obj._primer_check_forward(item1)
        self.assertEqual(item1, exp_item1)
        self.assertEqual(wf_obj.FinalState['Sequence'], 'CC')
        self.assertEqual(wf_obj.FinalState['Qual'], array([7,8]))
        self.assertEqual(wf_obj.FinalState['Forward primer'], 'AATTGG')
        self.assertFalse(wf_obj.Failed)

        wf_obj.wf_init(item2)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj.FinalState['Final barcode'] = 'AAAAAAAAAAAA'
        wf_obj._primer_check_forward(item2)
        self.assertEqual(item2, exp_item2)
        self.assertEqual(wf_obj.FinalState['Sequence'], 'CC')
        self.assertEqual(wf_obj.FinalState['Qual'], array([7,8]))
        self.assertEqual(wf_obj.FinalState['Forward primer'], 'AATTGC')
        self.assertFalse(wf_obj.Failed)

        wf_obj.wf_init(item3)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj.FinalState['Final barcode'] = 'AAAAAAAAAAAA'
        wf_obj._primer_check_forward(item3)
        self.assertEqual(item3, exp_item3)
        self.assertEqual(wf_obj.FinalState['Sequence'], 'GGTTGCCC')
        self.assertEqual(wf_obj.FinalState['Qual'], None)
        self.assertEqual(wf_obj.FinalState['Forward primer'], None)
        self.assertTrue(wf_obj.Failed)

        # item is not modified in place as retain priemr is True
        wf_obj = self._make_workflow_obj({'max_primer_mismatch':2,
                                          'retain_primer':True})
        item1 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'AATTGGCC',
                 'Qual':array([1,2,3,4,5,6,7,8])}
        item2 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'AATTGCCC',
                 'Qual':array([1,2,3,4,5,6,7,8])}
        item3 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'GGTTGCCC',
                 'Qual':array([1,2,3,4,5,6,7,8])}
        exp_item1 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'AATTGGCC',
                 'Qual':array([1,2,3,4,5,6,7,8])}
        exp_item2 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'AATTGCCC',
                 'Qual':array([1,2,3,4,5,6,7,8])}
        exp_item3 = {'Barcode':'AAAAAAAAAAAA', 'Sequence':'GGTTGCCC',
                 'Qual':array([1,2,3,4,5,6,7,8])}

        wf_obj.wf_init(item1)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj.FinalState['Final barcode'] = 'AAAAAAAAAAAA'
        wf_obj._primer_check_forward(item1)
        self.assertEqual(item1, exp_item1)
        self.assertEqual(wf_obj.FinalState['Sequence'], 'AATTGGCC')
        self.assertEqual(wf_obj.FinalState['Qual'], array([1,2,3,4,5,6,7,8]))
        self.assertEqual(wf_obj.FinalState['Forward primer'], 'AATTGG')
        self.assertFalse(wf_obj.Failed)

        wf_obj.wf_init(item2)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj.FinalState['Final barcode'] = 'AAAAAAAAAAAA'
        wf_obj._primer_check_forward(item2)
        self.assertEqual(item2, exp_item2)
        self.assertEqual(wf_obj.FinalState['Sequence'], 'AATTGCCC')
        self.assertEqual(wf_obj.FinalState['Qual'], array([1,2,3,4,5,6,7,8]))
        self.assertEqual(wf_obj.FinalState['Forward primer'], 'AATTGC')
        self.assertFalse(wf_obj.Failed)

        wf_obj.wf_init(item3)
        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj.FinalState['Final barcode'] = 'AAAAAAAAAAAA'
        wf_obj._primer_check_forward(item3)
        self.assertEqual(item3, exp_item3)
        self.assertEqual(wf_obj.FinalState['Sequence'], 'GGTTGCCC')
        self.assertEqual(wf_obj.FinalState['Qual'], None)
        self.assertEqual(wf_obj.FinalState['Forward primer'], None)
        self.assertTrue(wf_obj.Failed)

    def test_sequence_length_check(self):
        """Check the length of the sequence"""
        wf_obj = self._make_workflow_obj({'min_seq_len':5})
        item1 = {'Sequence':'AATTGGCC'}
        item2 = {'Sequence':'AATT'}

        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._sequence_length_check(item1)
        self.assertFalse(wf_obj.Failed)

        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._sequence_length_check(item2)
        self.assertTrue(wf_obj.Failed)

    def test_sequence_ambiguous_count(self):
        wf_obj = self._make_workflow_obj({'ambiguous_count':2})
        item1 = {'Sequence':'AATTGGCC'}
        item2 = {'Sequence':'AANNNTT'}
        item3 = {'Sequence':'AANTT'}

        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._sequence_ambiguous_count(item1)
        self.assertFalse(wf_obj.Failed)

        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._sequence_ambiguous_count(item2)
        self.assertTrue(wf_obj.Failed)

        wf_obj.Failed = False # note, normally handled by Workflow.__call__
        wf_obj._sequence_ambiguous_count(item3)
        self.assertFalse(wf_obj.Failed)

fasta1_simple = """>a
abcde
>b
asdasdasd
>c
123123
"""

fasta2_simple = """>x
abcdefg
>y
popopo
"""

qual1_simple = """>a
1 2 3 4 5
>b
1 1 1 1 1 1 1 1 1
>c
2 2 2 2 2 2
"""

qual2_simple = """>x
1 2 3 4 5 6 7
>y
1 1 1 1 1 1
"""

qual2_simple_bad = """>x
1 2 3 4 5 6
>y
1 1 1 1 1 1
"""

fastq1_simple = """@a
abcde
+a
abcde
@b
asdasdasd
+b
asdasdasd
@c
123123
+c
123123
"""

fastq2_simple = """@x
abcdefg
+x
abcdefg
@y
popopo
+y
popopo
"""

barcodes1_simple = """@a
test1
+a
1234
@b
test2
+b
12345
@c
test3
+c
aaccb
"""

barcodes2_simple = """@x
test4
+x
12312
@y
test5
+y
33333
"""

mapping = MetadataMap(
    {'s1':{'BarcodeSequence':'AAAAAAAAAAAA', 'LinkerPrimerSequence':'AATTGG,AATTCC'},
     's2':{'BarcodeSequence':'AAAAAAAAAAAC', 'LinkerPrimerSequence':''},
     's3':{'BarcodeSequence':'AAAAAAAAAAAG', 'LinkerPrimerSequence':''},
     's4':{'BarcodeSequence':'AAAAAAAAAAAT', 'LinkerPrimerSequence':''},
     's5':{'BarcodeSequence':'GGAGACAAGGGA', 'LinkerPrimerSequence':''}
    }, [])

fastq1 = """@990:2:4:11271:5323#1/1
GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC
+
bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`
@990:2:4:11271:5323#1/1
GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG
+
bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT
@990:2:4:11272:9538#1/1
GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG
+
b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa
@990:2:4:11272:9538#1/1
GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC
+
bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O
@990:2:4:11272:7447#1/1
GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGGCAGTCTAACCGCAAGGAGGACGCTGTCG
+
b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:7447#1/1
GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACAGCGTCCTCCTTGCTGTTAGACTTCCGGC
+
b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:19991#1/1
GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGGTGGGGCAACCGGCTGTCCCTTTTAGCGG
+
bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:19991#1/1
GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGGCTGGCCCTTTCCACCCA
+
bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOcBBBBBBBBBBBBBBBBB
@990:2:4:11272:4315#1/1
GTACTCACCGCCCGTCACGCCATGGGAGTTGGGCTTACCTGAAGCCCGCGAGCTAACCGGAAAGGGGGGGATGTGG
+
bbbb_bbbbbbbbbb```Q```BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:4315#1/1
GGCTACCTTGTTACGACTTCACCCCCGTCGCTCGGCGTACCTTCGACCGCTGCCTCCTTTTGGTTATATCTCCGGG
+
``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:5533#1/1
GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC
+
``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:5533#0/1
GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCGGCTGGCTCCCCTTTCGGGGGTACCTCAC
+
bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB
"""

barcode_fastq1 = """@990:2:4:11271:5323#1/2
AAAAAAAAAAAA
+
bbbbbbbbbbbb
@990:2:4:11271:5323#1/2
AAAAAAAAAAAC
+
bbcbbbbbbbbb
@990:2:4:11272:9538#1/2
AAAAAAAAAAAA
+
b_bbbbbbbbbb
@990:2:4:11272:9538#1/2
AAAAAAAAAAAT
+
bb^bbbbbbbbb
@990:2:4:11272:7447#1/2
AAAAAAAAAAAA
+
b`bbbbbbbbbb
@990:2:4:11272:7447#1/2
AAAAAAAAAAAA
+
b`bbbbbbbbbb
@990:2:4:11272:19991#1/2
AAAAAAAAAAAC
+
bbbbbbbbbbbb
@990:2:4:11272:19991#1/2
AAAAAAAAAAAC
+
bbbbbbbbbbbb
@990:2:4:11272:4315#1/2
AAAAAAAAAAAT
+
bbbb_bbbbbbb
@990:2:4:11272:4315#1/2
AAAAAAAAAAAT
+
``Q``````_``
@990:2:4:11272:5533#1/2
GAAAAAAAAAAT
+
``Q``````_``
@990:2:4:11272:5533#0/2
AAAAAAAAAAAT
+
bbbbbbbbbbbb
"""

fastq2 = """@M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0
GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC
+
bbbbbbbbbbBBBBBBBBBBBBBBBY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`
@M00176:17:000000000-A0CNA:1:1:17088:1773 1:N:0:0
GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG
+
bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT
@M00176:17:000000000-A0CNA:1:1:16738:1773 1:N:0:0
GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG
+
b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa
@M00176:17:000000000-A0CNA:1:1:12561:1773 1:N:0:0
GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC
+
bb^bbbBBBBbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O
@M00176:17:000000000-A0CNA:1:1:14596:1773 1:N:0:0
GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGGCAGTCTAACCGCAAGGAGGACGCTGTCG
+
b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`############################
@M00176:17:000000000-A0CNA:1:1:12515:1774 1:N:0:0
GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACAGCGTCCTCCTTGCTGTTAGACTTCCGGC
+
b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`############################
@M00176:17:000000000-A0CNA:1:1:17491:1774 1:N:0:0
GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGGTGGGGCAACCGGCTGTCCCTTTTAGCGG
+
bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb############################
@M00176:17:000000000-A0CNA:1:1:16427:1774 1:N:0:0
GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGGCTGGCCCTTTCCACCCA
+
bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc#################
@M00176:17:000000000-A0CNA:1:1:13372:1775 1:N:0:0
GTACTCACCGCCCGTCACGCCATGGGAGTTGGGCTTACCTGAAGCCCGCGAGCTAACCGGAAAGGGGGGGATGTGG
+
bbbb_bbbbbbbbbb```Q```######################################################
@M00176:17:000000000-A0CNA:1:1:14806:1775 1:N:0:0
GGCTACCTTGTTACGACTTCACCCCCGTCGCTCGGCGTACCTTCGACCGCTGCCTCCTTTTGGTTATATCTCCGGG
+
``Q``````_``BBBB````K]]a####################################################
@M00176:17:000000000-A0CNA:1:1:13533:1775 1:N:0:0
GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC
+
``Q``````_``````````K]]a####################################################
@M00176:17:000000000-A0CNA:1:1:18209:1775 1:N:0:0
GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCGGCTGGCTCCCCTTTCGGGGGTACCTCAC
+
bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb############################
"""

barcode_fastq2 = """@M00176:17:000000000-A0CNA:1:1:15487:1773 2:N:0:0
AAAAAAAAAAAA
+
bbbbbbbbbbbb
@M00176:17:000000000-A0CNA:1:1:17088:1773 2:N:0:0
AAAAAAAAAAAC
+
bbcbbbbbbbbb
@M00176:17:000000000-A0CNA:1:1:16738:1773 2:N:0:0
AAAAAAAAAAAA
+
b_bbbbbbbbbb
@M00176:17:000000000-A0CNA:1:1:12561:1773 2:N:0:0
AAAAAAAAAAAT
+
bb^bbbbbbbbb
@M00176:17:000000000-A0CNA:1:1:14596:1773 2:N:0:0
AAAAAAAAAAAA
+
b`bbbbbbbbbb
@M00176:17:000000000-A0CNA:1:1:12515:1774 2:N:0:0
AAAAAAAAAAAA
+
b`bbbbbbbbbb
@M00176:17:000000000-A0CNA:1:1:17491:1774 2:N:0:0
AAAAAAAAAAAC
+
bbbbbbbbbbbb
@M00176:17:000000000-A0CNA:1:1:16427:1774 2:N:0:0
AAAAAAAAAAAC
+
bbbbbbbbbbbb
@M00176:17:000000000-A0CNA:1:1:13372:1775 2:N:0:0
AAAAAAAAAAAT
+
bbbb_bbbbbbb
@M00176:17:000000000-A0CNA:1:1:14806:1775 2:N:0:0
AAAAAAAAAAAT
+
``Q``````_``
@M00176:17:000000000-A0CNA:1:1:13533:1775 2:N:0:0
GAAAAAAAAAAT
+
``Q``````_``
@M00176:17:000000000-A0CNA:1:1:18209:1775 2:N:0:0
AAAAAAAAAAAT
+
bbbbbbbbbbbb
"""


fastq1_expected_no_qual_unassigned = [
    ("s1_0 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`",
     0),
    ("s2_1 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     "bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT",
     1),
    ("s1_2 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     "b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa",
     2),
    ("s4_3 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     "bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O",
     3),
    ("s1_4 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGGCAGTCTAACCGCAAGGAGGACGCTGTCG",
     "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB",
     4),
    ("s1_5 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACAGCGTCCTCCTTGCTGTTAGACTTCCGGC",
     "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB",
     5),
    ("s2_6 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGGTGGGGCAACCGGCTGTCCCTTTTAGCGG",
     "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
     6),
    ("s2_7 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGGCTGGCCCTTTCCACCCA",
     "bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOcBBBBBBBBBBBBBBBBB",
     7),
    ("s4_8 990:2:4:11272:4315#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GTACTCACCGCCCGTCACGCCATGGGAGTTGGGCTTACCTGAAGCCCGCGAGCTAACCGGAAAGGGGGGGATGTGG",
     "bbbb_bbbbbbbbbb```Q```BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
     8),
    ("s4_9 990:2:4:11272:4315#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCCCGTCGCTCGGCGTACCTTCGACCGCTGCCTCCTTTTGGTTATATCTCCGGG",
     "``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
     9),
    ("Unassigned_10 990:2:4:11272:5533#1/1 orig_bc=GAAAAAAAAAAT new_bc=GAAAAAAAAAAT bc_diffs=0",
     "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
     "``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
     10),
    ("s4_11 990:2:4:11272:5533#0/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCGGCTGGCTCCCCTTTCGGGGGTACCTCAC",
     "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
     11)]

fastq1_expected_default = [
    ("s1_0 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`",
     0),
    ("s2_1 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     "bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT",
     1),
    ("s1_2 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     "b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa",
     2),
    ("s4_3 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     "bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O",
     3),
    ("s1_4 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGG",
     "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`",
     4),
    ("s1_5 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACA",
     "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`",
     5),
    ("s2_6 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGG",
     "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb",
     6),
    ("s2_7 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGG",
     "bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc",
     7),
    ("s4_8 990:2:4:11272:5533#0/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCG",
     "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb", 8)]

fastq1_expected_single_barcode = [
    ("s1_0 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`",
     0),
    ("s1_1 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     "bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT",
     1),
    ("s1_2 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     "b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa",
     2),
    ("s1_3 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     "bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O",
     3),
    ("s1_4 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGG",
     "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`",
     4),
    ("s1_5 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACA",
     "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`",
     5),
    ("s1_6 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGG",
     "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb",
     6),
    ("s1_7 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGG",
     "bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc",
     7),
    ("s1_8 990:2:4:11272:5533#0/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCG",
     "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb", 8)]

fastq2_expected_default = [
    ("s1_0 M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     "bbbbbbbbbbBBBBBBBBBBBBBBBY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`",
     0),
    ("s2_1 M00176:17:000000000-A0CNA:1:1:17088:1773 1:N:0:0 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     "bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT",
     1),
    ("s1_2 M00176:17:000000000-A0CNA:1:1:16738:1773 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     "b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa",
     2),
    ("s4_3 M00176:17:000000000-A0CNA:1:1:12561:1773 1:N:0:0 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     "bb^bbbBBBBbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O",
     3),
    ("s1_4 M00176:17:000000000-A0CNA:1:1:14596:1773 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGG",
     "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`",
     4),
    ("s1_5 M00176:17:000000000-A0CNA:1:1:12515:1774 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACA",
     "b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`",
     5),
    ("s2_6 M00176:17:000000000-A0CNA:1:1:17491:1774 1:N:0:0 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGG",
     "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb",
     6),
    ("s2_7 M00176:17:000000000-A0CNA:1:1:16427:1774 1:N:0:0 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGG",
     "bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc",
     7),
    ("s4_8 M00176:17:000000000-A0CNA:1:1:18209:1775 1:N:0:0 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCG",
     "bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb", 8)]

if __name__ == '__main__':
    main()
