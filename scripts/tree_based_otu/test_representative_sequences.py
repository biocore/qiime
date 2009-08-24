#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from representative_sequences import parse_otu_lines, load_tree,load_alignment,\
    RepresentativeSequence, LongestSequence, RandomSequence, \
    BestSharedSubstring, NoRepresentativeSequence, random_tiebreaker_f,\
    otus_to_seqids, RepSeqNode
from cogent.core.moltype import DNA
from cogent.core.alignment import Alignment
from cogent.parse.tree import DndParser

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Prototype"

class RepSeqNodeTests(TestCase):
    """Tests for RepSeqNode"""
    def test_setTipNames(self):
        """Sets TipNames attribute"""
        t_str = "((a,b)c,(d,e)f)g;"
        t = DndParser(t_str, constructor=RepSeqNode)
        t.setTipNames()
        self.assertEqualItems(t.TipNames, ['a','b','d','e'])
        self.assertEqualItems(t.Children[0].TipNames,['a','b'])
        self.assertEqualItems(t.Children[1].TipNames,['d','e'])
        self.assertEqualItems(t.getNodeMatchingName('a').TipNames, ['a'])
        self.assertEqualItems(t.getNodeMatchingName('b').TipNames, ['b'])
        self.assertEqualItems(t.getNodeMatchingName('d').TipNames, ['d'])
        self.assertEqualItems(t.getNodeMatchingName('e').TipNames, ['e'])

class RepresentativeSequenceTests(TestCase):
    """Tests for RepresentativeSequence baseclass"""
    def setUp(self):
        self.otu_map = {'1':['a','b'],
                        '2':['c'],
                        '3':['d','e','f']}
        self.seqs = {'a':'aaaatttt',
                     'b':'aaaatttg',
                     'c':'gggggggg',
                     'd':'ccccgggc',
                     'e':'cccc---c',
                     'f':'---cgggc'}
        self.aln = Alignment(self.seqs, MolType=DNA)

    def test_init(self):
        """Tests for RepresentativeSequence init"""
        rs = RepresentativeSequence(self.aln, self.otu_map)
        self.assertEqual(rs.Alignment, self.aln)
        self.assertEqual(rs.OTUmap, self.otu_map)
        self.assertEqual(rs.TieBreaker, random_tiebreaker_f)
        self.assertEqual(rs.Params, {})

    def test_call(self):
        """Tests for RepresentativeSequence call, should raise"""
        rs = RepresentativeSequence(self.aln, self.otu_map)
        self.assertRaises(NotImplementedError, rs) 

    def test_get_representative_seq_id(self):
        """Implemented in subclass, should raise"""
        rs = RepresentativeSequence(self.aln, self.otu_map)
        self.assertRaises(NotImplementedError, rs._get_representative_seq_id,\
                None) 

    def test_get_otu_cluster(self):
        """Returns OTU cluster, raises if none"""
        rs = RepresentativeSequence(self.aln, self.otu_map)
        exp = Alignment({'a':'aaaatttt','b':'aaaatttg'}, MolType=DNA)
        obs = rs._get_otu_cluster('1')

        rs.OTUmap['10'] = []
        self.assertRaises(NoRepresentativeSequence, rs._get_otu_cluster, '10')

class LongestSequenceTests(TestCase):
    """Tests for LongestSequence"""
    def setUp(self):
        self.otu_map = {'1':['a','b'],
                        '2':['c'],
                        '3':['d','e','f']}
        self.seqs = {'a':'aaaatttt',
                     'b':'aaaatttg',
                     'c':'gggggggg',
                     'd':'ccccgggc',
                     'e':'cccc---c',
                     'f':'---cgggc'}
        self.aln = Alignment(self.seqs, MolType=DNA)

    def test_get_representative_seq_id(self):
        """Should return expected sequence ids"""
        tiebreaker_f = lambda x,y: y[0]
        ls = LongestSequence(self.aln, self.otu_map, tiebreaker_f)

        exp = 'a'
        obs = ls._get_representative_seq_id('1')
        self.assertEqual(obs, exp)

        exp = 'c'
        obs = ls._get_representative_seq_id('2')
        self.assertEqual(obs, exp)

        exp = 'd'
        obs = ls._get_representative_seq_id('3')
        self.assertEqual(obs, exp)

    def test_call(self):
        """Should return SequenceCollection and map"""
        tiebreaker_f = lambda x,y: y[0]
        ls = LongestSequence(self.aln, self.otu_map, tiebreaker_f)
        exp_seqs = Alignment({'a':'aaaatttt','c':'gggggggg','d':'ccccgggc'})
        exp_map = {'1':'a','2':'c','3':'d'}

        obs_seqs, obs_map = ls()
        self.assertEqual(obs_seqs, exp_seqs)
        self.assertEqual(obs_map, exp_map)

class BestSharedSubstringTests(TestCase):
    def setUp(self):
        self.fail('not fully implemented')

class SupportMethodTests(TestCase):
    def test_parse_otu_lines(self):
        exp = {'0.01':['a','c','f'],
               '0.02':['a','q']}
        obs = parse_otu_lines(otu_lines)
        self.assertEqual(dict(obs), exp)

    def test_load_tree(self):
        t_str = "((a,b)c,(d,e)f)g;"
        obs_tree, obs_dict = load_tree(t_str)
        self.assertEqual(obs_tree.getNodeMatchingName('a').Id, '0')
        self.assertEqual(obs_tree.getNodeMatchingName('b').Id, '1')
        self.assertEqual(obs_tree.getNodeMatchingName('c').Id, '4')
        self.assertEqual(obs_tree.getNodeMatchingName('d').Id, '2')
        self.assertEqual(obs_tree.getNodeMatchingName('e').Id, '3')
        self.assertEqual(obs_tree.getNodeMatchingName('f').Id, '5')
        self.assertEqual(obs_tree.getNodeMatchingName('g').Id, '6')

        exp_map = {'0':obs_tree.getNodeMatchingName('a'),
                   '1':obs_tree.getNodeMatchingName('b'),
                   '2':obs_tree.getNodeMatchingName('d'),
                   '3':obs_tree.getNodeMatchingName('e'),
                   '4':obs_tree.getNodeMatchingName('c'),
                   '5':obs_tree.getNodeMatchingName('f'),
                   '6':obs_tree.getNodeMatchingName('g')}

        self.assertEqual(obs_dict, exp_map)

    def test_load_alignment(self):
        exp = Alignment({'a':'aaa--ttt','b':'agtttg--'})
        obs = load_alignment(aln_lines)
        self.assertEqual(obs, exp)

    def test_otus_to_seqids(self):
        """Returns mapping otus to descending seqids"""
        OTUs = [1,2]
        tree = DndParser("((a,b)c,(d,e)f)g;", constructor=RepSeqNode)
        node_dict = {1:tree.getNodeMatchingName('c'),
                     2:tree.getNodeMatchingName('f')}
        exp = {1:['a','b'],2:['d','e']}
        obs = otus_to_seqids(OTUs, node_dict)
        self.assertEqual(obs,exp)

otu_lines = """a\t0\t0.01
c\tNone\t0.01
f\t10\t0.01
a\t0\t0.02
q\t0\t0.02""".splitlines()

aln_lines = """>a foo bar
aaa.-ttt
>b bar foo
agtttg--""".splitlines()

if __name__ == '__main__':
    main()
