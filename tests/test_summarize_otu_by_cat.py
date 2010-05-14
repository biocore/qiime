#!/usr/bin/env python
#file test_summarize_otu_by_cat.py

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Julia Goodrich"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Julia Goodrich"
__email__ = "julia.goodrich@colorado.edu"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
from qiime.summarize_otu_by_cat import get_sample_cat_info, get_counts_by_cat

class TopLevelTests(TestCase):
    """Tests of top-level functions"""
    def setUp(self):
        """Set up for the tests in TopLevelTests Class"""
        self.map_file = """#SampleID\tDay\ttime\tDescription
#This is some comment about the study
1\t090809\t1200\tsome description of sample1
2\t090809\t1800\tsome description of sample2
3\t090909\t1200\tsome description of sample3
4\t090909\t1800\tsome description of sample4
5\t091009\t1200\tsome description of sample5"""
        self.cat_by_sample = {"1":[("Day","090809"),("time","1200"),('Description', 'some description of sample1')],
                              "2":[("Day","090809"),("time","1800"),('Description', 'some description of sample2')],
                              "3":[("Day","090909"),("time","1200"),('Description', 'some description of sample3')],
                              "4":[("Day","090909"),("time","1800"),('Description', 'some description of sample4')],
                              "5":[("Day","091009"),("time","1200"),('Description', 'some description of sample5')]}
        self.sample_by_cat = {("Day","090809"):["1","2"],
                              ("Day","090909"):["3","4"],
                              ("Day","091009"):["5"],
                              ("time","1200"):["1","3","5"],
                              ("time","1800"):["2","4"],
                              ('Description', 'some description of sample1'):['1'],
                              ('Description', 'some description of sample2'):['2'],
                              ('Description', 'some description of sample3'):['3'],
                              ('Description', 'some description of sample4'):['4'],
                              ('Description', 'some description of sample5'):['5']}
        self.num_cats = 3
        self.meta_dict = {"1":[("090809",0)],
                              "2":[("090809",0)],
                              "3":[("090909",0)],
                              "4":[("090909",0)],
                              "5":[("091009",0)]}
        self.labels_lists_dict = {"Day":["090809","090909","091009"],
                              "time":["1200","1800"],
                              'Description': ['some description of sample1', 'some description of sample2', 'some description of sample3', 'some description of sample4', 'some description of sample5']}
        self.num_samples_by_cat = {("Day","090809"):2,("Day","090909"):2,\
                                  ("Day","091009"):1,("time","1200"):3,\
                                  ("time","1800"):2,
                                    ('Description', 'some description of sample1'):1,
                                    ('Description', 'some description of sample2'):1,
                                    ('Description', 'some description of sample3'):1,
                                    ('Description', 'some description of sample4'):1,
                                    ('Description', 'some description of sample5'):1}
        self.otu_sample_file = """#Full OTU Counts
#OTU ID\t1\t2\t3\t4\t5\tConsensus Lineage
otu_1\t0\t1\t0\t0\t6\tBacteria; Actinobacteria; Coriobacteridae
otu_2\t2\t0\t0\t0\t0\tBacteria; Bacteroidetes; Bacteroidales; Bacteroidaceae
otu_3\t0\t0\t3\t1\t0\tBacteria; Firmicutes; Clostridia; Clostridiales
otu_4\t0\t0\t0\t0\t5\tBacteria; Spirochaetes; Spirochaetales; Spirochaetaceae
otu_5\t0\t4\t2\t0\t0\tBacteria; Bacteroidetes; Bacteroidales; Rikenellaceae
otu_6\t3\t6\t0\t0\t0\tBacteria; Bacteroidetes; Bacteroidales; Dysgonomonaceae
otu_7\t0\t0\t4\t2\t0\tBacteria; Bacteroidetes; Bacteroidales; Odoribacteriaceae
otu_8\t0\t0\t0\t0\t3\tBacteria; Bacteroidetes; Bacteroidales; Dysgonomonaceae; otu_425
otu_9\t2\t0\t0\t5\t0\tBacteria; Bacteroidetes; Bacteroidales; Dysgonomonaceae; otu_425
otu_10\t0\t2\t0\t4\t0\tBacteria; Firmicutes; Mollicutes; Clostridium_aff_innocuum_CM970"""
        self.cat_otu_table = """1\t0\t6
2\t0\t0
0\t4\t0
0\t0\t5
4\t2\t0
9\t0\t0
0\t6\t0
0\t0\t3
2\t5\t0
2\t4\t0"""


        self.cat_otu_table_norm = """3.84615\t0.0\t42.85714
14.28571\t0.0\t0.0
0.0\t20.83333\t0.0
0.0\t0.0\t35.71429
15.38462\t11.11111\t0.0
44.50549\t0.0\t0.0
0.0\t30.55556\t0.0
0.0\t0.0\t21.42857
14.28571\t20.83333\t0.0
7.69231\t16.66667\t0.0"""


        self.otus = ["otu_1","otu_2","otu_3","otu_4","otu_5","otu_6","otu_7","otu_8","otu_9","otu_10"]

        self.taxonomy = ["Bacteria; Actinobacteria; Coriobacteridae",
                      "Bacteria; Bacteroidetes; Bacteroidales; Bacteroidaceae",
                      "Bacteria; Firmicutes; Clostridia; Clostridiales",
                      "Bacteria; Spirochaetes; Spirochaetales; Spirochaetaceae",
                      "Bacteria; Bacteroidetes; Bacteroidales; Rikenellaceae",
                      "Bacteria; Bacteroidetes; Bacteroidales; Dysgonomonaceae",
                      "Bacteria; Bacteroidetes; Bacteroidales; Odoribacteriaceae",
                      "Bacteria; Bacteroidetes; Bacteroidales; Dysgonomonaceae; otu_425",
                      "Bacteria; Bacteroidetes; Bacteroidales; Dysgonomonaceae; otu_425",
                      "Bacteria; Firmicutes; Mollicutes; Clostridium_aff_innocuum_CM970"]
        self.num_cat = {"Day":2,
                        "time":4}

    def test_get_sample_cat_info(self):
        """get_sample_cat_info should return hand calculated values"""
        cat_by_sample, sample_by_cat, num_cats, meta_dict,label_lists_dict,num_samples_by_cat = get_sample_cat_info(self.map_file.split('\n'),"Day")
        self.assertEqual(cat_by_sample,self.cat_by_sample)
        self.assertEqual(sample_by_cat,self.sample_by_cat)
        self.assertEqual(num_cats,self.num_cats)
        self.assertEqual(meta_dict,self.meta_dict)
        self.assertEqual(label_lists_dict,self.labels_lists_dict)
        self.assertEqual(num_samples_by_cat,self.num_samples_by_cat)

    def test_get_counts_by_cat(self):
        """get_counts_by_cat should return hand calculated values"""
        cat_otu_table, otus, taxonomy = get_counts_by_cat(\
               self.otu_sample_file.split('\n'), self.num_cats,\
                self.meta_dict,self.labels_lists_dict["Day"],"Day",self.num_samples_by_cat,False)
        cat_otu_table_test = []
        for l in cat_otu_table:
            cat_otu_table_test.append('\t'.join(map(str,l)))
        self.assertEqual('\n'.join(cat_otu_table_test),self.cat_otu_table)
        self.assertEqual(otus,self.otus)
        self.assertEqual(taxonomy,self.taxonomy)
        cat_otu_table, otus, taxonomy = get_counts_by_cat(\
                self.otu_sample_file.split('\n'), self.num_cats,\
                self.meta_dict,self.labels_lists_dict["Day"],"Day",self.num_samples_by_cat,True)
        cat_otu_table_test = []
        for l in cat_otu_table:
            cat_otu_table_test.append('\t'.join(map(str,l)))
        self.assertEqual('\n'.join(cat_otu_table_test),self.cat_otu_table_norm)
        self.assertEqual(otus,self.otus)
        self.assertEqual(taxonomy,self.taxonomy)

if __name__ =='__main__':
    main()
