#!/usr/bin/env python
#file test_split_libraries.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

from cogent.util.unit_test import TestCase, main
from StringIO import StringIO
from numpy import array
from pipe454.split_libraries import (
    generate_possibilities, get_infile, count_mismatches,
    ok_mm_primer, check_map, get_primer_seqs, fasta_ids, qual_score,
    qual_scores, count_ambig, split_seq, primer_exceeds_mismatches,
    check_barcode, make_histograms, format_histograms, SeqQualBad
)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_generate_possibilities(self):
        """generate_possibilities should make possible strings"""
        self.assertEqual(generate_possibilities('ACG'), ['ACG'])
        self.assertEqual(generate_possibilities('RGY'), 
            ['AGT', 'AGC', 'GGT', 'GGC'])

    def test_get_infile(self):
        """get_infile should return filehandle"""
        pass    #not practically testable, but obvious file I/O

    def test_count_mismatches(self):
        """count_mismatches should count mismatches correctly"""
        self.assertEqual(count_mismatches('GG','GG',10), 0)
        self.assertEqual(count_mismatches('GGG','AAA',10), 3)
        self.assertEqual(count_mismatches('GGG','AAA',1), 2)

    def test_ok_mm_primer(self):
        """ok_mm_primer should test that primer is within max mismatches"""
        primers = ['AAAA', 'GGGG']
        self.assertEqual(ok_mm_primer('AAAA', primers, 0), True)
        self.assertEqual(ok_mm_primer('AAAA', primers, 3), True)
        self.assertEqual(ok_mm_primer('CCCC', primers, 0), False)
        self.assertEqual(ok_mm_primer('CCCA', primers, 3), True)
        self.assertEqual(ok_mm_primer('CCCA', primers, 2), False)
        self.assertEqual(ok_mm_primer('CCGG', primers, 2), True)
        self.assertEqual(ok_mm_primer('CCGA', primers, 2), False)

    def test_check_map(self):
        """check_map should return valid barcodes as expected"""
        s = """#SampleID\tBarcodeSequence\tX\tDescription
#fake data
x\tAA\t3\tsample x
y\t"AC"\t4\t"sample y"
z\tGG\t5\tsample z"""
        f = StringIO(s)
        f.name='test.xls'
        headers, id_map, barcode_to_sample_id, warnings, errors = \
            check_map(f)

        self.assertEqual(barcode_to_sample_id, {'AA':'x','AC':'y','GG':'z'})
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])

    def test_get_primer_seqs(self):
        """get_primer_seqs should resolve degen patterns, check length"""
        primer_pats = ['AGY', 'GAG']
        res = get_primer_seqs(primer_pats, min_primer_len=2)
        self.assertEqual(res, (set(['AGT','AGC','GAG']),3))
        self.assertRaises(ValueError, get_primer_seqs, primer_pats, 5)
        self.assertRaises(ValueError, get_primer_seqs, primer_pats+['x'], 3)

    def test_fasta_ids(self):
        """fasta_ids should return list of ids in fasta files, no dups"""
        first = StringIO('>x\nACT\n>y\nAAA')
        first_copy = StringIO('>x\nACT\n>y\nAAA')
        second = StringIO('>a\nGGG\n>b\nCCC')
        self.assertEqual(fasta_ids([first, second]), set(['x','y','a','b']))
        first.seek(0) #need to reset so we can read it again
        self.assertRaises(ValueError, fasta_ids, [first,first_copy])

    def test_qual_score(self):
        """qual_score should return dict of {id: qual_scores}"""
        scores = StringIO('>x\n5 10 5\n12\n>y\n30 40')
        self.assertEqual(qual_score(scores), {'x':[5,10,5,12],'y':[30,40]})

    def test_qual_scores(self):
        """qual_scores should return dict of {id:qual_scores}"""
        scores = StringIO('>x\n5 10 5\n12\n>y\n30 40')
        scores2= StringIO('>a\n5 10 5\n12\n>b\n30 40')
        self.assertEqual(qual_scores([scores, scores2]),
            {'x':[5,10,5,12],'y':[30,40],'a':[5,10,5,12],'b':[30,40]})

    def test_count_ambig(self):
        """count_ambig should count ambiguous bases in seq"""
        s = 'ACC'
        s2 = 'RNY'
        s3 = 'NA'
        self.assertEqual(count_ambig(s), 0)
        self.assertEqual(count_ambig(s2), 3)
        self.assertEqual(count_ambig(s3), 1)
        self.assertEqual(count_ambig(''), 0)

    def test_split_seq(self):
        """split_seq should split seq into pieces"""
        seq = 'AAAACCCCCGTGTGTGT'
        barcode, primer, remainder = split_seq(seq, 4, 5)
        self.assertEqual(barcode, 'AAAA')
        self.assertEqual(primer, 'CCCCC')
        self.assertEqual(remainder, 'GTGTGTGT')

    def test_primer_exceeds_mismatches(self):
        """primer_exceeds_mismatches returns True if too many mismatches"""
        primers = ['AAAA', 'TTTT']
        exact = 'AAAA'
        mismatch_ok = 'AAAT'
        mismatch_bad = 'GGGG'
        self.assertEqual(primer_exceeds_mismatches(exact, primers, 0), False)
        self.assertEqual(primer_exceeds_mismatches(mismatch_ok, primers, 1), 
            False)
        self.assertEqual(primer_exceeds_mismatches(mismatch_bad, primers, 2), 
            True)

    def test_check_barcode(self):
        """check_barcode should return False if barcode ok, True otherwise"""
        self.assertEqual(check_barcode('AA', None, ['AA']), (False, 'AA'))
        self.assertEqual(check_barcode('GCATCGTCCACA', 'golay_12', 
            ['GCATCGTCAACA']), (True, 'GCATCGTCAACA'))
        self.assertEqual(check_barcode('AA', None, ['TT']), (True, None))

    def test_make_histograms(self):
        """make_histograms should make correct histograms"""
        pre_lengths = [100, 110, 105, 130, 135]
        post_lengths = [130, 135, 150]
        pre_hist, post_hist, bin_edges = \
            make_histograms(pre_lengths, post_lengths)
        self.assertEqual(pre_hist, array([2,1,0,2,0,0]))
        self.assertEqual(post_hist, array([0,0,0,2,0,1]))
        self.assertEqual(bin_edges, array([100,110,120,130,140,150,160]))

    def test_format_histograms(self):
        """format_histograms should print histograms correctly"""
        self.assertEqual(format_histograms(array([2,1,0,2,0,0]),
            array([0,0,0,2,0,1]), array([100,110,120,130,140,150,160])),
            """Length\tBefore\tAfter
100\t2\t0
110\t1\t0
120\t0\t0
130\t2\t2
140\t0\t0
150\t0\t1""")

class SeqQualBadTests(TestCase):
    """Tests of the SeqQualBad class"""
    def test_init(self):
        """SeqQualBad should init OK"""
        sq = SeqQualBad('Q', None)
        self.assertEqual(sq.Name, 'Q')
        self.assertEqual(sq.F, None)
        self.assertEqual(sq.FailedIds, [])

    def test_call(self):
        """SeqQualBad should apply correct function"""
        f = lambda id_, seq, qual: len(seq) > 3
        s1 = 'aa'
        s2 = 'aaaa'
        sq = SeqQualBad('Q', f)
        self.assertEqual(sq('x',s1, [1,2,3]), False)
        self.assertEqual(sq('y',s2, [1,2,3]), True)
        self.assertEqual(sq.FailedIds, ['y'])

    def test_str(self):
        """SeqQualBad should apply correct function"""
        f = lambda id_, seq, qual: len(seq) > 3
        s1 = 'aa'
        s2 = 'aaaa'
        sq = SeqQualBad('Q', f)
        self.assertEqual(sq('x',s1, [1,2,3]), False)
        self.assertEqual(str(sq), 'Q\t0')
        self.assertEqual(sq('y',s2, [1,2,3]), True)
        self.assertEqual(str(sq), 'Q\t1')


if __name__ =='__main__':
    main()
