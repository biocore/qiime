#!/usr/bin/env python

from itertools import chain, izip
from numpy import array
from cogent.util.unit_test import TestCase, main
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.fastq import MinimalFastqParser
from qiime.parse import MinimalQualParser
from qiime.process_seqs import (_fasta_qual_gen,
        fasta_iterator, _fastq_barcode_gen, fastq_iterator, _fastq_gen,
        SequenceWorkflow)
from qiime.quality import ascii_to_phred64
from qiime.util import MetadataMap

class FastqIteratorTests(TestCase):
    def setUp(self):
        fastq1_gen = MinimalFastqParser(fastq1_simple.splitlines())
        fastq2_gen = MinimalFastqParser(fastq2_simple.splitlines())
        barcodes1_gen = MinimalFastqParser(barcodes1_simple.splitlines())
        barcodes2_gen = MinimalFastqParser(barcodes2_simple.splitlines())
        
        self.fastq_gen = chain(fastq1_gen, fastq2_gen)
        self.barcodes_gen = chain(barcodes1_gen, barcodes2_gen)

        self.reversed_fastq_gen = chain(fastq2_gen, fastq1_gen)

    def test_fastq_barcode_gen_simple(self):
        exp_data = [('a', 'abcde', 'test1', array([33,34,35,36,37])),
                    ('b', 'asdasdasd', 'test2', array([33,51,36] * 3)),
                    ('c', '123123', 'test3', array([-15, -14, -13] * 2)),
                    ('x', 'abcdefg', 'test4', array([33,34,35,36,37,38,39])),
                    ('y', 'popopo', 'test5', array([48,47] * 3))]
        exp = []
        for id_,seq,bc,qual in exp_data:
            exp.append({'SequenceID':id_, 'Sequence':seq, 'Qual':qual, 
                        'Barcode':bc})
        
        obs = _fastq_barcode_gen(self.fastq_gen, self.barcodes_gen,
                                 ascii_to_phred64)
        for o,e in izip(obs,exp):
            self.assertEqual(o['SequenceID'], e['SequenceID'])
            self.assertEqual(o['Sequence'], e['Sequence'])
            self.assertTrue((o['Qual'] == e['Qual']).all())
            self.assertEqual(o['Barcode'], e['Barcode'])

    def test_fasta_barcode_gen_mismatch_ids(self):
        with self.assertRaises(ValueError):
            g = _fasta_qual_gen(self.reversed_fastq_gen, self.barcodes_gen)
            _ = list(g)
    
    def test_fastq_iterators_just_fastq(self):
        exp_data = [('a', 'abcde', array([33,34,35,36,37])),
                    ('b', 'asdasdasd', array([33,51,36] * 3)),
                    ('c', '123123', array([-15, -14, -13] * 2)),
                    ('x', 'abcdefg', array([33,34,35,36,37,38,39])),
                    ('y', 'popopo', array([48,47] * 3))]
        exp = []
        for id_,seq,qual in exp_data:
            exp.append({'SequenceID':id_, 'Sequence':seq, 'Qual':qual, 
                        'Barcode':None})
        
        open_fps = map(lambda x: x.splitlines(), [fastq1_simple, fastq2_simple])
        obs = [d.copy() for d in fastq_iterator(open_fps)]
        self.assertEqual(obs, exp)

    def test_fastq_iterators_barcodes(self):
        exp_data = [('a', 'abcde', 'test1', array([33,34,35,36,37])),
                    ('b', 'asdasdasd', 'test2', array([33,51,36] * 3)),
                    ('c', '123123', 'test3', array([-15, -14, -13] * 2)),
                    ('x', 'abcdefg', 'test4', array([33,34,35,36,37,38,39])),
                    ('y', 'popopo', 'test5', array([48,47] * 3))]
        exp = []
        for id_,seq,bc,qual in exp_data:
            exp.append({'SequenceID':id_, 'Sequence':seq, 'Qual':qual, 
                        'Barcode':bc})
        
        splitter = lambda x: x.splitlines()
        fastq_fps = map(splitter, [fastq1_simple, fastq2_simple])
        bc_fps = map(splitter, [barcodes1_simple, barcodes2_simple])

        obs = fastq_iterator(fastq_fps, bc_fps)
        for o,e in izip(obs,exp):
            self.assertEqual(o['SequenceID'], e['SequenceID'])
            self.assertEqual(o['Sequence'], e['Sequence'])
            self.assertTrue((o['Qual'] == e['Qual']).all())
            self.assertEqual(o['Barcode'], e['Barcode'])
        
class FastaIteratorTests(TestCase):
    def setUp(self):
        fasta1_gen = MinimalFastaParser(fasta1_simple.splitlines())
        qual1_gen = MinimalQualParser(qual1_simple.splitlines())
        fasta2_gen = MinimalFastaParser(fasta2_simple.splitlines())
        qual2_gen = MinimalQualParser(qual2_simple.splitlines())
        qual2_bad_gen = MinimalQualParser(qual2_simple_bad.splitlines())

        self.fasta_gen = chain(fasta1_gen, fasta2_gen)
        self.qual_gen = chain(qual1_gen, qual2_gen)

        self.reversed_fasta_gen = chain(fasta2_gen, fasta1_gen)
        self.qual_bad_gen = chain(qual1_gen, qual2_bad_gen)

    def test_fasta_qual_gen_simple(self):
        exp_data = [('a', 'abcde', array([1, 2, 3, 4, 5])),
                    ('b', 'asdasdasd', array([1,1,1,1,1,1,1,1,1])),
                    ('c', '123123', array([2, 2, 2, 2, 2, 2])),
                    ('x', 'abcdefg', array([1, 2, 3, 4, 5, 6, 7])),
                    ('y', 'popopo', array([1, 1, 1, 1, 1, 1]))]
        exp = []
        for id_,seq,qual in exp_data:
            exp.append({'SequenceID':id_, 'Sequence':seq, 'Qual':qual, 
                        'Barcode':None})
        
        obs = _fasta_qual_gen(self.fasta_gen, self.qual_gen)
        for o,e in izip(obs,exp):
            self.assertEqual(o['SequenceID'], e['SequenceID'])
            self.assertEqual(o['Sequence'], e['Sequence'])
            self.assertTrue((o['Qual'] == e['Qual']).all())

    def test_fasta_qual_gen_mismatch_ids(self):
        with self.assertRaises(ValueError):
            g = _fasta_qual_gen(self.reversed_fasta_gen, self.qual_gen)
            _ = list(g)
    
    def test_fasta_qual_gen_mismatch_length(self):
        with self.assertRaises(ValueError):
            _ = list(_fasta_qual_gen(self.fasta_gen, self.qual_bad_gen))

    def test_fasta_iterators_just_fasta(self):
        exp_data = [('a', 'abcde', None),
                    ('b', 'asdasdasd', None),
                    ('c', '123123', None),
                    ('x', 'abcdefg', None),
                    ('y', 'popopo', None)]

        exp = []
        for id_,seq,qual in exp_data:
            exp.append({'SequenceID':id_, 'Sequence':seq, 'Qual':qual, 
                        'Barcode':None})
        
        open_fps = map(lambda x: x.splitlines(), [fasta1_simple, fasta2_simple])
        obs = [d.copy() for d in fasta_iterator(open_fps)]
        self.assertEqual(obs, exp)

    def test_fasta_iterators_fasta_qual(self):
        exp_data = [('a', 'abcde', array([1, 2, 3, 4, 5])),
                    ('b', 'asdasdasd', array([1,1,1,1,1,1,1,1,1])),
                    ('c', '123123', array([2, 2, 2, 2, 2, 2])),
                    ('x', 'abcdefg', array([1, 2, 3, 4, 5, 6, 7])),
                    ('y', 'popopo', array([1, 1, 1, 1, 1, 1]))]

        exp = []
        for id_,seq,qual in exp_data:
            exp.append({'SequenceID':id_, 'Sequence':seq, 'Qual':qual, 
                        'Barcode':None})
        splitter = lambda x: x.splitlines()
        fasta_fps = map(splitter, [fasta1_simple, fasta2_simple])
        qual_fps = map(splitter, [qual1_simple, qual2_simple])

        obs = fasta_iterator(fasta_fps, qual_fps)
        for o,e in izip(obs, exp):
            self.assertEqual(o['SequenceID'], e['SequenceID'])
            self.assertEqual(o['Sequence'], e['Sequence'])
            self.assertTrue((o['Qual'] == e['Qual']).all())


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
        return SequenceWorkflow(Options=options, Mapping=self.mapping)

    def test_workflow_construction(self):
        x = self._make_workflow_obj({'foo':'bar'})
    
    def test_wf_init(self):
        wf_obj = self._make_workflow_obj({'foo':'bar'})
        wf_obj.FinalState['Sequence'] = 'w00t'
        wf_obj.wf_init({'Sequence':'foo'})
        self.assertEqual(set(wf_obj.FinalState.values()), set([None, 'foo']))

    def test_quality_max_bad_run_length(self):
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

        wf_obj._quality_min_per_read_length_fraction(item1)
        self.assertFalse(wf_obj.Failed)

        wf_obj._quality_min_per_read_length_fraction(item2)
        self.assertFalse(wf_obj.Failed)

        wf_obj._quality_min_per_read_length_fraction(item3)
        self.assertTrue(wf_obj.Failed)

        self.assertEqual(item1, exp1)
        self.assertEqual(item2, exp2)
        self.assertEqual(item3, exp3)
        
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
    {'s1':{'BarcodeSequence':'AAAAAAAAAAAA', 'LinkerPrimerSequence':''},
     's2':{'BarcodeSequence':'AAAAAAAAAAAC', 'LinkerPrimerSequence':''},
     's3':{'BarcodeSequence':'AAAAAAAAAAAG', 'LinkerPrimerSequence':''},
     's4':{'BarcodeSequence':'AAAAAAAAAAAT', 'LinkerPrimerSequence':''}
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
