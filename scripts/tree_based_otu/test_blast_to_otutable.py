#!/usr/bin/env python

from numpy import array
from cogent.util.unit_test import TestCase, main
from blast_to_otutable import convert_to_matrix, get_best_hit, clean_queryid,\
        best_by_evalue, parse_seq_to_db, parse_seq_to_otu

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Prototype"

class TestCases(TestCase):
    def test_convert_to_matrix(self):
        input = {('s1','o1'):5,
                 ('s1','o2'):1,
                 ('s2','o3'):10,
                 ('s2','o1'):3}
        exp_matrix = array([[5,3],[1,0],[0,10]])
        exp_sampleids = ['s1','s2']
        exp_otus = ['o1','o2','o3']
        
        obs_matrix, obs_otus, obs_sampleids = convert_to_matrix(input)
        
        self.assertEqual(obs_matrix, exp_matrix)
        self.assertEqual(obs_otus, exp_otus)
        self.assertEqual(obs_sampleids, exp_sampleids)

    def test_get_best_hit(self):
        record = ('>asdasd', [{'a':None,'E-VALUE':0.5,'SUBJECT ID':10},
                              {'b':'TACO','E-VALUE':0.2,'SUBJECT ID':15},
                              {'c':1,'E-VALUE':0.6,'SUBJECT ID':7}])
        exp = ('asdasd',15)
        obs = get_best_hit(record)
        self.assertEqual(obs, exp)

    def test_clean_queryid(self):
        foo = '>asdasd'
        exp = 'asdasd'
        obs = clean_queryid(foo)
        self.assertEqual(obs, exp)

    def test_best_by_evalue(self):
        recs = [{'a':None,'E-VALUE':0.5,'SUBJECT ID':10},
                {'b':'TACO','E-VALUE':0.2,'SUBJECT ID':15},
                {'c':1,'E-VALUE':0.6,'SUBJECT ID':7}]
        exp = (0.2, {'b':'TACO','E-VALUE':0.2,'SUBJECT ID':15})
        obs = best_by_evalue(recs)
        self.assertEqual(obs, exp)

    def test_parse_seq_to_db(self):
        exp = {'1':'a_b','2':'d_e','3':'g_h'}
        obs = parse_seq_to_db(seq_to_db_mapping_lines)
        self.assertEqual(obs, exp)

    def test_parse_seq_to_otu(self):
        """Parses Sequence id to otu mapping file"""
        exp = {'1':'2','3':'4','5':'6'}
        obs = parse_seq_to_otu(seq_to_otu_mapping_lines)
        self.assertEqual(obs, exp)

seq_to_otu_mapping_lines = """1\t2
3\t4
5\t6""".splitlines()

seq_to_db_mapping_lines = """1\ta\tb\tc
2\td\te\tf
3\tg\th\ti""".splitlines()

if __name__ == '__main__':
    main()
