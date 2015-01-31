#!/usr/bin/env python
# file test_make_otu_network.py

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2011, The QIIME Project"  # consider project name
# remember to add yourself
__credits__ = ["Julia Goodrich", "Jose Carlos Clemente Litran"]
__license__ = "GPL"
__version__ = "1.9.0-dev"
__maintainer__ = "Jose Clemente"
__email__ = "jose.clemente@gmail.com"

from os import remove, close
from os.path import exists
import os
import shutil
from tempfile import mkstemp

from numpy import array
from random import choice, randrange
from tempfile import mkdtemp
from unittest import TestCase, main
from qiime.stats import G_2_by_2
from qiime.make_otu_network import get_sample_info, get_connection_info, \
    get_num_con_cat, get_num_cat, make_table_file, make_stats_files,\
    make_props_files
from qiime.util import load_qiime_config, write_biom_table
from biom.table import Table


class OtuNetworkTests(TestCase):

    def setUp(self):
        self.qiime_config = load_qiime_config()
        self.tmp_dir = self.qiime_config['temp_dir'] or '/tmp/'

        self.map_file = """#SampleID	Day	time	Description
#This is some comment about the study
1	090809	1200	some description of sample1
2	090809	1800	some description of sample2
3	090909	1200	some description of sample3
4	090909	1800	some description of sample4
5	091009	1200	some description of sample5"""
        self.cat_by_sample = {"1": [("Day", "090809"), ("time", "1200")],
                              "2": [("Day", "090809"), ("time", "1800")],
                              "3": [("Day", "090909"), ("time", "1200")],
                              "4": [("Day", "090909"), ("time", "1800")],
                              "5": [("Day", "091009"), ("time", "1200")]}
        self.sample_by_cat = {("Day", "090809"): ["1", "2"],
                              ("Day", "090909"): ["3", "4"],
                              ("Day", "091009"): ["5"],
                              ("time", "1200"): ["1", "3", "5"],
                              ("time", "1800"): ["2", "4"]}

        self.num_cats = 2
        self.meta_dict = {"1": ["090809	1200", 0],
                          "2": ["090809	1800", 0],
                          "3": ["090909	1200", 0],
                          "4": ["090909	1800", 0],
                          "5": ["091009	1200", 0]}
        self.labels = ["from", "to", "eweight", "consensus_lin", "Day", "time"]
        self.node_labels = ["node_name", "node_disp_name", "ntype", "degree",
                            "weighted_degree", "consensus_lin", "Day", "time"]
        self.label_list = [["090809", "090909", "091009"], ["1200", "1800"]]

        self.otu_table_vals = array([[0, 1, 0, 0, 6],
                                     [2, 0, 0, 0, 0],
                                     [0, 0, 3, 1, 0],
                                     [0, 0, 0, 0, 5],
                                     [0, 4, 2, 0, 0],
                                     [3, 6, 0, 0, 0],
                                     [0, 0, 4, 2, 0],
                                     [0, 0, 0, 0, 3],
                                     [2, 0, 0, 5, 0],
                                     [0, 2, 0, 4, 0]])

        otu_table = Table(self.otu_table_vals,
                          ['otu_1', 'otu_2', 'otu_3', 'otu_4', 'otu_5',
                           'otu_6', 'otu_7', 'otu_8', 'otu_9', 'otu_10'],
                          ['1', '2', '3', '4', '5'],
                          [{"taxonomy": ["Bacteria", "Actinobacteria",
                                         "Coriobacteridae"]},
                           {"taxonomy": ["Bacteria", "Bacteroidetes",
                                         "Bacteroidales", "Bacteroidaceae"]},
                           {"taxonomy": ["Bacteria", "Firmicutes",
                                         "Clostridia", "Clostridiales"]},
                           {"taxonomy": ["Bacteria", "Spirochaetes",
                                         "Spirochaetales", "Spirochaetaceae"]},
                           {"taxonomy": ["Bacteria", "Bacteroidetes",
                                         "Bacteroidales", "Rikenellaceae"]},
                           {"taxonomy": ["Bacteria", "Bacteroidetes",
                                         "Bacteroidales", "Dysgonomonaceae"]},
                           {"taxonomy": ["Bacteria", "Bacteroidetes",
                                         "Bacteroidales",
                                         "Odoribacteriaceae"]},
                           {"taxonomy": ["Bacteria", "Bacteroidetes",
                                         "Bacteroidales", "Dysgonomonaceae",
                                         "otu_425"]},
                           {"taxonomy": ["Bacteria", "Bacteroidetes",
                                         "Bacteroidales", "Dysgonomonaceae",
                                         "otu_425"]},
                           {"taxonomy": ["Bacteria", "Firmicutes",
                                         "Mollicutes",
                                         "Clostridium_aff_innocuum_CM970"]}],
                          [None, None, None, None, None])

        fd, self.otu_table_fp = mkstemp(
            dir=self.tmp_dir, prefix='test_make_otu_network_otu_table',
            suffix='.biom')
        close(fd)
        write_biom_table(otu_table, self.otu_table_fp)

        self.otu_sample_file = """#Full OTU Counts
#OTU ID	1	2	3	4	5	Consensus Lineage
otu_1	0	1	0	0	6	Bacteria; Actinobacteria; Coriobacteridae
otu_2	2	0	0	0	0	Bacteria; Bacteroidetes; Bacteroidales; Bacteroidaceae
otu_3	0	0	3	1	0	Bacteria; Firmicutes; Clostridia; Clostridiales
otu_4	0	0	0	0	5	Bacteria; Spirochaetes; Spirochaetales; Spirochaetaceae
otu_5	0	4	2	0	0	Bacteria; Bacteroidetes; Bacteroidales; Rikenellaceae
otu_6	3	6	0	0	0	Bacteria; Bacteroidetes; Bacteroidales; Dysgonomonaceae
otu_7	0	0	4	2	0	Bacteria; Bacteroidetes; Bacteroidales; Odoribacteriaceae
otu_8	0	0	0	0	3	Bacteria; Bacteroidetes; Bacteroidales; Dysgonomonaceae; otu_425
otu_9	2	0	0	5	0	Bacteria; Bacteroidetes; Bacteroidales; Dysgonomonaceae; otu_425
otu_10	0	2	0	4	0	Bacteria; Firmicutes; Mollicutes; Clostridium_aff_innocuum_CM970"""

        self.con_by_sample = {
            '1': set(['2', '4']), '2': set(['5', '3', '1', '4']),
            '3': set(['4', '2']), '4': set(['3', '1', '2']),
            '5': set(['2'])}

        self.edge_file_str = [
            "2	otu_1	1.0	Bacteria:Actinobacteria:Coriobacteridae	090809	1800",
            "5	otu_1	6.0	Bacteria:Actinobacteria:Coriobacteridae	091009	1200",
            "1	otu_2	2.0	Bacteria:Bacteroidetes:Bacteroidales:Bacteroidaceae	090809	1200",
            "3	otu_3	3.0	Bacteria:Firmicutes:Clostridia:Clostridiales	090909	1200",
            "4	otu_3	1.0	Bacteria:Firmicutes:Clostridia:Clostridiales	090909	1800",
            "5	otu_4	5.0	Bacteria:Spirochaetes:Spirochaetales:Spirochaetaceae	091009	1200",
            "2	otu_5	4.0	Bacteria:Bacteroidetes:Bacteroidales:Rikenellaceae	090809	1800",
            "3	otu_5	2.0	Bacteria:Bacteroidetes:Bacteroidales:Rikenellaceae	090909	1200",
            "1	otu_6	3.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae	090809	1200",
            "2	otu_6	6.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae	090809	1800",
            "3	otu_7	4.0	Bacteria:Bacteroidetes:Bacteroidales:Odoribacteriaceae	090909	1200",
            "4	otu_7	2.0	Bacteria:Bacteroidetes:Bacteroidales:Odoribacteriaceae	090909	1800",
            "5	otu_8	3.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae:otu_425	091009	1200",
            "1	otu_9	2.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae:otu_425	090809	1200",
            "4	otu_9	5.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae:otu_425	090909	1800",
            "2	otu_10	2.0	Bacteria:Firmicutes:Mollicutes:Clostridium_aff_innocuum_CM970	090809	1800",
            "4	otu_10	4.0	Bacteria:Firmicutes:Mollicutes:Clostridium_aff_innocuum_CM970	090909	1800"]

        self.node_file_str = ["1	1	user_node	3	7.0	other	090809	1200",
                              "2	2	user_node	4	13.0	other	090809	1800",
                              "3	3	user_node	3	9.0	other	090909	1200",
                              "4	4	user_node	4	12.0	other	090909	1800",
                              "5	5	user_node	3	14.0	other	091009	1200",
                              "otu_1		otu_node	2	7.0	Bacteria:Actinobacteria:Coriobacteridae	otu	otu",
                              "otu_2		otu_node	1	2.0	Bacteria:Bacteroidetes:Bacteroidales:Bacteroidaceae	otu	otu",
                              "otu_3		otu_node	2	4.0	Bacteria:Firmicutes:Clostridia:Clostridiales	otu	otu",
                              "otu_4		otu_node	1	5.0	Bacteria:Spirochaetes:Spirochaetales:Spirochaetaceae	otu	otu",
                              "otu_5		otu_node	2	6.0	Bacteria:Bacteroidetes:Bacteroidales:Rikenellaceae	otu	otu",
                              "otu_6		otu_node	2	9.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae	otu	otu",
                              "otu_7		otu_node	2	6.0	Bacteria:Bacteroidetes:Bacteroidales:Odoribacteriaceae	otu	otu",
                              "otu_8		otu_node	1	3.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae:otu_425	otu	otu",
                              "otu_9		otu_node	2	7.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae:otu_425	otu	otu",
                              "otu_10		otu_node	2	6.0	Bacteria:Firmicutes:Mollicutes:Clostridium_aff_innocuum_CM970	otu	otu"]

        self.red_edge_file_str = [
            "2	otu_1	1.0	Bacteria:Actinobacteria:Coriobacteridae	090809	1800",
            "5	otu_1	6.0	Bacteria:Actinobacteria:Coriobacteridae	091009	1200",
            "1	@1	1.0	missed	090809	1200",
            "3	otu_3	3.0	Bacteria:Firmicutes:Clostridia:Clostridiales	090909	1200",
            "4	otu_3	1.0	Bacteria:Firmicutes:Clostridia:Clostridiales	090909	1800",
            "5	@5	1.0	missed	091009	1200",
            "2	otu_5	4.0	Bacteria:Bacteroidetes:Bacteroidales:Rikenellaceae	090809	1800",
            "3	otu_5	2.0	Bacteria:Bacteroidetes:Bacteroidales:Rikenellaceae	090909	1200",
            "1	otu_6	3.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae	090809	1200",
            "2	otu_6	6.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae	090809	1800",
            "3	otu_7	4.0	Bacteria:Bacteroidetes:Bacteroidales:Odoribacteriaceae	090909	1200",
            "4	otu_7	2.0	Bacteria:Bacteroidetes:Bacteroidales:Odoribacteriaceae	090909	1800",
            "1	otu_9	2.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae:otu_425	090809	1200",
            "4	otu_9	5.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae:otu_425	090909	1800",
            "2	otu_10	2.0	Bacteria:Firmicutes:Mollicutes:Clostridium_aff_innocuum_CM970	090809	1800",
            "4	otu_10	4.0	Bacteria:Firmicutes:Mollicutes:Clostridium_aff_innocuum_CM970	090909	1800"]

        self.red_node_file_str = ["1	1	user_node	3	7.0	other	090809	1200",
                                  "2	2	user_node	4	13.0	other	090809	1800",
                                  "3	3	user_node	3	9.0	other	090909	1200",
                                  "4	4	user_node	4	12.0	other	090909	1800",
                                  "5	5	user_node	3	14.0	other	091009	1200",
                                  "otu_1		otu_node	2	7.0	Bacteria:Actinobacteria:Coriobacteridae	otu	otu",
                                  "@1		otu_collapsed	1	1.0	other	otu	otu",
                                  "otu_3		otu_node	2	4.0	Bacteria:Firmicutes:Clostridia:Clostridiales	otu	otu",
                                  "@5		otu_collapsed	2	2.0	other	otu	otu",
                                  "otu_5		otu_node	2	6.0	Bacteria:Bacteroidetes:Bacteroidales:Rikenellaceae	otu	otu",
                                  "otu_6		otu_node	2	9.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae	otu	otu",
                                  "otu_7		otu_node	2	6.0	Bacteria:Bacteroidetes:Bacteroidales:Odoribacteriaceae	otu	otu",
                                  "otu_9		otu_node	2	7.0	Bacteria:Bacteroidetes:Bacteroidales:Dysgonomonaceae:otu_425	otu	otu",
                                  "otu_10		otu_node	2	6.0	Bacteria:Firmicutes:Mollicutes:Clostridium_aff_innocuum_CM970	otu	otu"]

        self.otu_dc = {1: 3, 2: 7}
        self.sample_dc = {3: 3, 4: 2}
        self.degree_counts = {1: 3, 2: 7, 3: 3, 4: 2}

        self.num_con_cat = {"Day": 2,
                            "time": 1}
        self.num_con = 6
        self.num_cat = {"Day": 2,
                        "time": 4}
        self.num_cat_less = {"Day": 1,
                             "time": 3}
        self._paths_to_clean_up = [self.otu_table_fp]
        self._dir_to_clean_up = ''

    def tearDown(self):
        map(remove, self._paths_to_clean_up)
        if self._dir_to_clean_up != '':
            shutil.rmtree(self._dir_to_clean_up)

    def test_get_sample_info(self):
        cat_by_sample, sample_by_cat, num_cats, meta_dict, labels, node_labels, \
            label_list = get_sample_info(self.map_file.split('\n'))
        self.assertEqual(cat_by_sample, self.cat_by_sample)
        self.assertEqual(sample_by_cat, self.sample_by_cat)
        self.assertEqual(num_cats, self.num_cats)
        self.assertEqual(meta_dict, self.meta_dict)
        self.assertEqual(labels, self.labels)
        self.assertEqual(node_labels, self.node_labels)
        self.assertEqual(label_list, self.label_list)

    def test_get_connection_info(self):
        con_by_sample, node_file_str, edge_file_str, red_node_file_str,\
            red_edge_file_str, otu_dc, degree_counts, sample_dc = \
            get_connection_info(self.otu_table_fp, self.num_cats,
                                self.meta_dict)
#           get_connection_info(self.otu_sample_file.split('\n'), self.num_cats,\
#                            self.meta_dict)

        self.assertEqual(con_by_sample, self.con_by_sample)
        self.assertEqual(set(node_file_str), set(self.node_file_str))
        self.assertEqual(set(edge_file_str), set(self.edge_file_str))
        self.assertEqual(set(red_node_file_str), set(self.red_node_file_str))
        self.assertEqual(set(red_edge_file_str), set(self.red_edge_file_str))
        self.assertEqual(otu_dc, self.otu_dc)
        self.assertEqual(degree_counts, self.degree_counts)
        self.assertEqual(sample_dc, self.sample_dc)

    def test_get_num_con_cat(self):
        num_con_cat, num_con = get_num_con_cat(
            self.con_by_sample, self.cat_by_sample)
        self.assertEqual(num_con_cat, self.num_con_cat)
        self.assertEqual(num_con, self.num_con)

    def test_get_num_cat(self):
        num_cat = get_num_cat(self.sample_by_cat, self.cat_by_sample.keys())
        self.assertEqual(num_cat, self.num_cat)
        num_cat = get_num_cat(
            self.sample_by_cat,
            self.cat_by_sample.keys()[:-1])
        self.assertEqual(num_cat, self.num_cat_less)

    def test_make_table_file(self):
        random_dir_name = mkdtemp()
        foldername = random_dir_name

        self._dir_to_clean_up = foldername

        try:
            os.mkdir(foldername)
        except OSError:
            pass

        obs = foldername

        try:
            os.mkdir(os.path.join(obs, "otu_network"))
        except OSError:
            pass

        obs = os.path.join(obs, "otu_network")
        make_table_file(
            self.edge_file_str,
            self.labels,
            obs,
            "real_edge_table.txt")

        self.assertTrue(exists(foldername + "/otu_network/real_edge_table.txt"), 'The file was not created in \
the appropriate location')

    def test_make_stats_files(self):
        random_dir_name = mkdtemp()
        foldername = random_dir_name
        self._dir_to_clean_up = foldername

        try:
            os.mkdir(foldername)
        except OSError:
            pass

        obs = foldername

        try:
            os.mkdir(os.path.join(obs, "otu_network"))
        except OSError:
            pass

        try:
            os.mkdir(os.path.join(obs, "otu_network/stats"))
        except OSError:
            pass

        obs = os.path.join(obs, "otu_network")
        make_stats_files(
            self.sample_dc,
            self.otu_dc,
            self.degree_counts,
            self.num_con_cat,
            self.num_con,
            self.num_cat,
            self.cat_by_sample,
            obs)

        self.assertTrue(exists(foldername + "/otu_network/stats/real_dc_otu_degree.txt"), 'The file was not created in \
the appropriate location')
        self.assertTrue(exists(foldername + "/otu_network/stats/real_dc_sample_degree.txt"), 'The file was not created in \
the appropriate location')
        self.assertTrue(exists(foldername + "/otu_network/stats/real_dc_sample_otu_degree.txt"), 'The file was not created in \
the appropriate location')
        self.assertTrue(exists(foldername + "/otu_network/stats/real_cat_stats_Day.txt"), 'The file was not created in \
the appropriate location')

        self.assertTrue(exists(foldername + "/otu_network/stats/real_cat_stats_time.txt"), 'The file was not created in \
the appropriate location')

    def test_make_props_files(self):
        random_dir_name = mkdtemp()
        foldername = random_dir_name

        self._dir_to_clean_up = foldername

        try:
            os.mkdir(foldername)
        except OSError:
            pass

        obs = foldername

        try:
            os.mkdir(os.path.join(obs, "otu_network"))
        except OSError:
            pass

        try:
            os.mkdir(os.path.join(obs, "otu_network/props"))
        except OSError:
            pass

        obs = os.path.join(obs, "otu_network")

        self.assertTrue(exists(foldername + "/otu_network/props/"), 'The file was not created in \
the appropriate location')
        self.assertTrue(exists(foldername + "/otu_network/props"), 'The file was not created in \
the appropriate location')


if __name__ == '__main__':
    main()
