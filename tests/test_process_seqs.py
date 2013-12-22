#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.fastq import MinimalFastqParser
from qiime.parse import MinimalQualParser
from itertools import chain
from numpy import array
from qiime.process_seqs import (_fasta_qual_gen,
        fasta_iterator, _fastq_barcode_gen, fastq_iterator, _fastq_gen)
from qiime.quality import ascii_to_phred64
from itertools import izip

class FastqIteratorTests(TestCase):
    def setUp(self):
        fastq1_gen = MinimalFastqParser(fastq1.splitlines())
        fastq2_gen = MinimalFastqParser(fastq2.splitlines())
        barcodes1_gen = MinimalFastqParser(barcodes1.splitlines())
        barcodes2_gen = MinimalFastqParser(barcodes2.splitlines())
        
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
        
        open_fps = map(lambda x: x.splitlines(), [fastq1, fastq2])
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
        fastq_fps = map(splitter, [fastq1, fastq2])
        bc_fps = map(splitter, [barcodes1, barcodes2])

        obs = fastq_iterator(fastq_fps, bc_fps)
        for o,e in izip(obs,exp):
            self.assertEqual(o['SequenceID'], e['SequenceID'])
            self.assertEqual(o['Sequence'], e['Sequence'])
            self.assertTrue((o['Qual'] == e['Qual']).all())
            self.assertEqual(o['Barcode'], e['Barcode'])
        
class FastaIteratorTests(TestCase):
    def setUp(self):
        fasta1_gen = MinimalFastaParser(fasta1.splitlines())
        qual1_gen = MinimalQualParser(qual1.splitlines())
        fasta2_gen = MinimalFastaParser(fasta2.splitlines())
        qual2_gen = MinimalQualParser(qual2.splitlines())
        qual2_bad_gen = MinimalQualParser(qual2_bad.splitlines())

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
        
        open_fps = map(lambda x: x.splitlines(), [fasta1, fasta2])
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
        fasta_fps = map(splitter, [fasta1, fasta2])
        qual_fps = map(splitter, [qual1, qual2])

        obs = fasta_iterator(fasta_fps, qual_fps)
        for o,e in izip(obs, exp):
            self.assertEqual(o['SequenceID'], e['SequenceID'])
            self.assertEqual(o['Sequence'], e['Sequence'])
            self.assertTrue((o['Qual'] == e['Qual']).all())

fasta1 = """>a
abcde
>b
asdasdasd
>c
123123
"""

fasta2 = """>x
abcdefg
>y
popopo
"""

qual1 = """>a
1 2 3 4 5
>b
1 1 1 1 1 1 1 1 1
>c
2 2 2 2 2 2
"""

qual2 = """>x
1 2 3 4 5 6 7
>y
1 1 1 1 1 1
"""

qual2_bad = """>x
1 2 3 4 5 6
>y
1 1 1 1 1 1
"""

fastq1 = """@a
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

fastq2 = """@x
abcdefg
+x
abcdefg
@y
popopo
+y
popopo
"""

barcodes1 = """@a
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

barcodes2 = """@x
test4
+x
12312
@y
test5
+y
33333
"""

if __name__ == '__main__':
    main()
