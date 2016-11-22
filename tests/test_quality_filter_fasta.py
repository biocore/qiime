#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import MinimalQualParser
from itertools import chain
from numpy import array
from qiime.quality_filter_fasta import (_fasta_qual_strict,
        fasta_qual_iterator)

class IteratorTests(TestCase):
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

    def test_fasta_qual_strict_simple(self):
        exp = [('a', 'abcde', 'a', array([1, 2, 3, 4, 5])),
               ('b', 'asdasdasd', 'b', array([1,1,1,1,1,1,1,1,1])),
               ('c', '123123', 'c', array([2, 2, 2, 2, 2, 2])),
               ('x', 'abcdefg', 'x', array([1, 2, 3, 4, 5, 6, 7])),
               ('y', 'popopo', 'y', array([1, 1, 1, 1, 1, 1]))]

        obs = _fasta_qual_strict(self.fasta_gen, self.qual_gen)
        for o,e in zip(obs,exp):
            osi, osd, oqi, oqd = o
            esi, esd, eqi, eqd = e
            self.assertEqual((osi, osd, oqi), (esi, esd, eqi))
            self.assertTrue((oqd == eqd).all())
    
    def test_fasta_qual_strict_mismatch_ids(self):
        with self.assertRaises(ValueError):
            g = _fasta_qual_strict(self.reversed_fasta_gen, self.qual_gen)
            _ = list(g)
    
    def test_fasta_qual_strict_mismatch_length(self):
        with self.assertRaises(ValueError):
            _ = list(_fasta_qual_strict(self.fasta_gen, self.qual_bad_gen))

    def test_fasta_qual_iterators_just_fasta(self):
        exp = [('a', 'abcde', None, None),
               ('b', 'asdasdasd', None, None),
               ('c', '123123', None, None),
               ('x', 'abcdefg', None, None),
               ('y', 'popopo', None, None)]

        open_fps = map(lambda x: x.splitlines(), [fasta1, fasta2])
        obs = list(fasta_qual_iterator(open_fps))
        self.assertEqual(obs, exp)

    def test_fasta_qual_iterators_fasta_qual(self):
        exp = [('a', 'abcde', 'a', array([1, 2, 3, 4, 5])),
               ('b', 'asdasdasd', 'b', array([1,1,1,1,1,1,1,1,1])),
               ('c', '123123', 'c', array([2, 2, 2, 2, 2, 2])),
               ('x', 'abcdefg', 'x', array([1, 2, 3, 4, 5, 6, 7])),
               ('y', 'popopo', 'y', array([1, 1, 1, 1, 1, 1]))]

        splitter = lambda x: x.splitlines()
        fasta_fps = map(splitter, [fasta1, fasta2])
        qual_fps = map(splitter, [qual1, qual2])

        obs = fasta_qual_iterator(fasta_fps, qual_fps)
        for o,e in zip(obs, exp):
            osi, osd, oqi, oqd = o
            esi, esd, eqi, eqd = e
            self.assertEqual((osi, osd, oqi), (esi, esd, eqi))
            self.assertTrue((oqd == eqd).all())

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

if __name__ == '__main__':
    main()
