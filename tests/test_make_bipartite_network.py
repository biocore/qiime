#!/usr/bin/env python
# File created on 17 Apr 2013
from __future__ import division

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"


from shutil import rmtree
from os.path import exists, join
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files

from qiime.util import get_qiime_temp_dir
from qiime.test import initiate_timeout, disable_timeout
from qiime.make_bipartite_network import (make_otu_node_table, make_edge_table,
                                          make_sample_node_table, make_node_attr_table)
from biom.parse import parse_biom_table
from qiime.parse import parse_mapping_file_to_dict


class TopLevelTests(TestCase):

    def setUp(self):
        '''Nothing required by all scripts. Pass.'''
        pass

    def test_make_sample_node_table(self):
        '''Test that the sample node table is created correctly.'''
        # test when sampleids in biom == sampleids in mapping file
        bt = parse_biom_table(BIOM_STRING_1)
        mf_dict = parse_mapping_file_to_dict(MF_LINES.split('\n'))[0]
        obs = make_sample_node_table(bt, mf_dict)
        exp = ['#NodeID\tNodeType\tAbundance\tTimePt\tStudy\tTreatment\tDiet',
               's1\tsample\t148.0\t1\ta\tpre\thf',
               's2\tsample\t156.0\t2\ta\tpre\tlf',
               's3\tsample\t164.0\t3\ta\tpre\thf',
               's4\tsample\t172.0\t4\ta\tpost\tlf',
               's5\tsample\t180.0\t5\ta\tpost\tmf']
        self.assertEqual(obs, exp)
        # test when sampleids in biom are a subset of sampleids in mapping file
        bt = parse_biom_table(BIOM_STRING_2)
        obs = make_sample_node_table(bt, mf_dict)
        exp = ['#NodeID\tNodeType\tAbundance\tTimePt\tStudy\tTreatment\tDiet',
               's3\tsample\t164.0\t3\ta\tpre\thf',
               's4\tsample\t172.0\t4\ta\tpost\tlf',
               's5\tsample\t180.0\t5\ta\tpost\tmf']
        self.assertEqual(obs, exp)

    def test_make_otu_node_table(self):
        '''Test that make_otu_node_table makes accurate calculations.'''
        # test when length of md_fields and length of split taxonomy wouldn't
        # agree when md type is list or string.
        bt1 = parse_biom_table(BIOM_STRING_1)
        md_key = 'taxonomy'
        md_fields = ['k']
        obs = make_otu_node_table(bt1, md_key, md_fields)
        exp = [
            '#NodeID\tNodeType\tAbundance\tk',
            'o1\totu\t15.0\tk__Bacteria',
            'o2\totu\t40.0\tk__Bacteria',
            'o3\totu\t65.0\tk__Bacteria',
            'o4\totu\t90.0\tk__Bacteria',
            'o5\totu\t115.0\tk__Bacteria',
            'o6\totu\t140.0\tk__Bacteria',
            'o7\totu\t165.0\tk__Bacteria',
            'o8\totu\t190.0\tk__Bacteria']
        self.assertEqual(obs, exp)
        md_fields = ['k', '1', '2', '3', '4', '5', '6']
        obs = make_otu_node_table(bt1, md_key, md_fields)
        exp = [
            '#NodeID\tNodeType\tAbundance\tk\t1\t2\t3\t4\t5\t6',
            'o1\totu\t15.0\tk__Bacteria\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae\tOther\tOther',
            'o2\totu\t40.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1\tOther\tOther',
            'o3\totu\t65.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1\tOther\tOther',
            'o4\totu\t90.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2\tOther\tOther',
            'o5\totu\t115.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2\tOther\tOther',
            'o6\totu\t140.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae\tOther\tOther',
            'o7\totu\t165.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae\tOther\tOther',
            'o8\totu\t190.0\tk__Bacteria\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae\tOther\tOther']

        # test when the length of the md_fields is correct
        md_fields = ['k', 'p', 'c', 'o', 'f']
        obs_bt1 = make_otu_node_table(bt1, md_key, md_fields)
        exp_bt1 = \
            ['#NodeID\tNodeType\tAbundance\tk\tp\tc\to\tf',
             'o1\totu\t15.0\tk__Bacteria\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae',
             'o2\totu\t40.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1',
             'o3\totu\t65.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1',
             'o4\totu\t90.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2',
             'o5\totu\t115.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2',
             'o6\totu\t140.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae',
             'o7\totu\t165.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae',
             'o8\totu\t190.0\tk__Bacteria\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae']
        self.assertEqual(obs_bt1, exp_bt1)
        md_fields = ['k', 'p', 'c', 'o', 'f']
        bt2 = parse_biom_table(BIOM_STRING_3)
        obs_bt2 = make_otu_node_table(bt2, md_key, md_fields)
        exp_bt2 = \
            ['#NodeID\tNodeType\tAbundance\tk\tp\tc\to\tf',
             'o1\totu\t12.0\tk__Bacteria\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae',
             'o2\totu\t27.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1',
             'o3\totu\t42.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1',
             'o4\totu\t57.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2',
             'o5\totu\t72.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2',
             'o6\totu\t87.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae',
             'o7\totu\t102.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae',
             'o8\totu\t117.0\tk__Bacteria\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae']
        self.assertEqual(obs_bt2, exp_bt2)

        # test when the md is of type dict and fields are correct
        bt = parse_biom_table(BIOM_STRING_4)
        md_fields = ['kingdom', 'phylum', 'class', 'order', 'family']
        obs_bt = make_otu_node_table(bt, md_key, md_fields)
        exp_bt = \
            ['#NodeID\tNodeType\tAbundance\tkingdom\tphylum\tclass\torder\tfamily',
             'o1\totu\t15.0\tk__Bacteria\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae',
             'o2\totu\t40.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1',
             'o3\totu\t65.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1',
             'o4\totu\t90.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2',
             'o5\totu\t115.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2',
             'o6\totu\t140.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae',
             'o7\totu\t165.0\tk__Bacteria\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae',
             'o8\totu\t190.0\tk__Bacteria\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae']
        self.assertEqual(obs_bt, exp_bt)
        # test that it raises an error when md_fields not found in md_dict
        md_fields = ['KINGDOM', 'phylum', 'class', 'order', 'family']
        self.assertRaises(
            ValueError,
            make_otu_node_table,
            bt,
            md_key,
            md_fields)
        # test that it doesn't error when md_fields is subset of md_dict
        md_fields = ['phylum', 'class', 'order', 'family']
        obs_bt = make_otu_node_table(bt, md_key, md_fields)
        exp_bt = \
            ['#NodeID\tNodeType\tAbundance\tphylum\tclass\torder\tfamily',
             'o1\totu\t15.0\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae',
             'o2\totu\t40.0\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1',
             'o3\totu\t65.0\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1',
             'o4\totu\t90.0\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2',
             'o5\totu\t115.0\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2',
             'o6\totu\t140.0\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae',
             'o7\totu\t165.0\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae',
             'o8\totu\t190.0\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae']
        self.assertEqual(obs_bt, exp_bt)

        # test when the md is type defaultdict
        bt = parse_biom_table(BIOM_STRING_5)
        # test that it raises an error when md_fields not found in md_dict
        md_fields = ['KINGDOM', 'phylum', 'class', 'order', 'family']
        self.assertRaises(
            ValueError,
            make_otu_node_table,
            bt,
            md_key,
            md_fields)
        # test that it doesn't error when md_fields is subset of md_dict
        md_fields = ['phylum', 'class', 'order', 'family']
        obs_bt = make_otu_node_table(bt, md_key, md_fields)
        exp_bt = \
            ['#NodeID\tNodeType\tAbundance\tphylum\tclass\torder\tfamily',
             'o1\totu\t15.0\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae',
             'o2\totu\t40.0\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1',
             'o3\totu\t65.0\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos1',
             'o4\totu\t90.0\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2',
             'o5\totu\t115.0\tp__Firmicutes\tc__Clostridia\to__Clostridiales\tf__Lachnos2',
             'o6\totu\t140.0\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae',
             'o7\totu\t165.0\tp__Firmicutes\tc__Clostridia\to__Clostri3\tf__Lachnospiraceae',
             'o8\totu\t190.0\tp__Firmicutes\tc__Clost5\to__Clostridiales\tf__Lachnospiraceae']
        self.assertEqual(obs_bt, exp_bt)

    def test_node_attr_table(self):
        '''Test that node_attr_table is made correctly.'''
        # test color,size,shape with shared attributes: NodeType, Abundance
        bt = parse_biom_table(BIOM_STRING_1)
        mf_dict = parse_mapping_file_to_dict(MF_LINES.split('\n'))[0]
        sample_node_lines = make_sample_node_table(bt, mf_dict)
        md_key = 'taxonomy'
        md_fields = ['k', 'p', 'c', 'o', 'f']
        otu_node_lines = make_otu_node_table(bt, md_key, md_fields)
        scolor = ['NodeType']
        ocolor = ['NodeType']
        sshape = ['NodeType', 'Abundance']
        oshape = ['NodeType']
        ssize = ['Abundance']
        osize = ['Abundance', 'NodeType']
        obs_out = make_node_attr_table(otu_node_lines, sample_node_lines,
                                       scolor, ocolor, ssize, osize, sshape, oshape)
        exp_out = \
            ['#NodeID\tNodeType\tAbundance\tColor\tSize\tShape',
             's1\tsample\t148.0\tsample\t148.0\tsample_148.0',
             's2\tsample\t156.0\tsample\t156.0\tsample_156.0',
             's3\tsample\t164.0\tsample\t164.0\tsample_164.0',
             's4\tsample\t172.0\tsample\t172.0\tsample_172.0',
             's5\tsample\t180.0\tsample\t180.0\tsample_180.0',
             'o1\totu\t15.0\totu\t15.0_otu\totu',
             'o2\totu\t40.0\totu\t40.0_otu\totu',
             'o3\totu\t65.0\totu\t65.0_otu\totu',
             'o4\totu\t90.0\totu\t90.0_otu\totu',
             'o5\totu\t115.0\totu\t115.0_otu\totu',
             'o6\totu\t140.0\totu\t140.0_otu\totu',
             'o7\totu\t165.0\totu\t165.0_otu\totu',
             'o8\totu\t190.0\totu\t190.0_otu\totu']
        # order different because computed by hand. As long as sets the same
        # we are confident they are the same
        self.assertEqual(set(obs_out), set(exp_out))

        # test color,size,shape with some shared and some non-shared attrs
        md_key = 'taxonomy'
        md_fields = ['k', 'p', 'c', 'o', 'f']
        otu_node_lines = make_otu_node_table(bt, md_key, md_fields)
        scolor = ['NodeType', 'Diet']
        ocolor = ['k', 'p']
        sshape = ['Treatment', 'Abundance']
        oshape = ['NodeType', 'o']
        ssize = ['Abundance']
        osize = ['Abundance', 'NodeType']

        obs_out = make_node_attr_table(otu_node_lines, sample_node_lines,
                                       scolor, ocolor, ssize, osize, sshape, oshape)
        exp_out = \
            ['#NodeID\tNodeType\tAbundance\tColor\tSize\tShape',
             's1\tsample\t148.0\tsample_hf\t148.0\tpre_148.0',
             's2\tsample\t156.0\tsample_lf\t156.0\tpre_156.0',
             's3\tsample\t164.0\tsample_hf\t164.0\tpre_164.0',
             's4\tsample\t172.0\tsample_lf\t172.0\tpost_172.0',
             's5\tsample\t180.0\tsample_mf\t180.0\tpost_180.0',
             'o1\totu\t15.0\tk__Bacteria_p__Firmicutes\t15.0_otu\totu_o__Clostridiales',
             'o2\totu\t40.0\tk__Bacteria_p__Firmicutes\t40.0_otu\totu_o__Clostridiales',
             'o3\totu\t65.0\tk__Bacteria_p__Firmicutes\t65.0_otu\totu_o__Clostridiales',
             'o4\totu\t90.0\tk__Bacteria_p__Firmicutes\t90.0_otu\totu_o__Clostridiales',
             'o5\totu\t115.0\tk__Bacteria_p__Firmicutes\t115.0_otu\totu_o__Clostridiales',
             'o6\totu\t140.0\tk__Bacteria_p__Firmicutes\t140.0_otu\totu_o__Clostri3',
             'o7\totu\t165.0\tk__Bacteria_p__Firmicutes\t165.0_otu\totu_o__Clostri3',
             'o8\totu\t190.0\tk__Bacteria_p__Firmicutes\t190.0_otu\totu_o__Clostridiales']
        # order different because computed by hand. As long as sets the same
        # we are confident they are the same
        self.assertEqual(set(obs_out), set(exp_out))

    def test_make_edge_table(self):
        '''Test that edge table is created properly.'''
        bt = parse_biom_table(BIOM_STRING_3)
        obs_out = make_edge_table(bt)
        exp_out = [
            '#Sample\tOTU\tAbundance',
            's3\to1\t3.0',
            's3\to2\t8.0',
            's3\to3\t13.0',
            's3\to4\t18.0',
            's3\to5\t23.0',
            's3\to6\t28.0',
            's3\to7\t33.0',
            's3\to8\t38.0',
            's4\to1\t4.0',
            's4\to2\t9.0',
            's4\to3\t14.0',
            's4\to4\t19.0',
            's4\to5\t24.0',
            's4\to6\t29.0',
            's4\to7\t34.0',
            's4\to8\t39.0',
            's5\to1\t5.0',
            's5\to2\t10.0',
            's5\to3\t15.0',
            's5\to4\t20.0',
            's5\to5\t25.0',
            's5\to6\t30.0',
            's5\to7\t35.0',
            's5\to8\t40.0']
        self.assertEqual(set(obs_out), set(exp_out))

        # test with a row and a col that are all zero
        bt = parse_biom_table(BIOM_STRING_6)
        obs_out = make_edge_table(bt)
        exp_out = [
            '#Sample\tOTU\tAbundance',
            's2\to1\t2.0',
            's2\to2\t7.0',
            's2\to3\t12.0',
            's2\to4\t17.0',
            's2\to6\t27.0',
            's2\to7\t32.0',
            's2\to8\t37.0',
            's3\to1\t3.0',
            's3\to2\t8.0',
            's3\to3\t13.0',
            's3\to4\t18.0',
            's3\to6\t28.0',
            's3\to7\t33.0',
            's3\to8\t38.0',
            's4\to1\t4.0',
            's4\to2\t9.0',
            's4\to3\t14.0',
            's4\to4\t19.0',
            's4\to6\t29.0',
            's4\to7\t34.0',
            's4\to8\t39.0',
            's5\to1\t5.0',
            's5\to2\t10.0',
            's5\to3\t15.0',
            's5\to4\t20.0',
            's5\to6\t30.0',
            's5\to7\t35.0',
            's5\to8\t40.0']
        self.assertEqual(set(obs_out), set(exp_out))

BIOM_STRING_1 = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "will","date": "2013-04-17T16:45:57.531405","matrix_type": "sparse","matrix_element_type": "float","shape": [8, 5],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,4.0],[0,4,5.0],[1,0,6.0],[1,1,7.0],[1,2,8.0],[1,3,9.0],[1,4,10.0],[2,0,11.0],[2,1,12.0],[2,2,13.0],[2,3,14.0],[2,4,15.0],[3,0,16.0],[3,1,17.0],[3,2,18.0],[3,3,19.0],[3,4,20.0],[4,0,21.0],[4,1,22.0],[4,2,23.0],[4,3,24.0],[4,4,25.0],[5,0,26.0],[5,1,27.0],[5,2,28.0],[5,3,29.0],[5,4,30.0],[6,0,31.0],[6,1,32.0],[6,2,33.0],[6,3,34.0],[6,4,35.0],[7,0,36.0],[7,1,37.0],[7,2,38.0],[7,3,39.0],[7,4,40.0]],"rows": [{"id": "o1", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clost5", "o__Clostridiales", "f__Lachnospiraceae"]}},{"id": "o2", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnos1"]}},{"id": "o3", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnos1"]}},{"id": "o4", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnos2"]}},{"id": "o5", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnos2"]}},{"id": "o6", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostri3", "f__Lachnospiraceae"]}},{"id": "o7", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostri3", "f__Lachnospiraceae"]}},{"id": "o8", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clost5", "o__Clostridiales", "f__Lachnospiraceae"]}}],"columns": [{"id": "s1", "metadata": null},{"id": "s2", "metadata": null},{"id": "s3", "metadata": null},{"id": "s4", "metadata": null},{"id": "s5", "metadata": null}]}'
# BIOM_STRING_2 lacks samples 1 and 2
BIOM_STRING_2 = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "will","date": "2013-04-17T19:04:21.047065","matrix_type": "sparse","matrix_element_type": "float","shape": [8, 3],"data": [[0,0,3.0],[0,1,4.0],[0,2,5.0],[1,0,8.0],[1,1,9.0],[1,2,10.0],[2,0,13.0],[2,1,14.0],[2,2,15.0],[3,0,18.0],[3,1,19.0],[3,2,20.0],[4,0,23.0],[4,1,24.0],[4,2,25.0],[5,0,28.0],[5,1,29.0],[5,2,30.0],[6,0,33.0],[6,1,34.0],[6,2,35.0],[7,0,38.0],[7,1,39.0],[7,2,40.0]],"rows": [{"id": "o1", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clost5", "o__Clostridiales", "f__Lachnospiraceae"]}},{"id": "o2", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnos1"]}},{"id": "o3", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnos1"]}},{"id": "o4", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnos2"]}},{"id": "o5", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnos2"]}},{"id": "o6", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostri3", "f__Lachnospiraceae"]}},{"id": "o7", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostri3", "f__Lachnospiraceae"]}},{"id": "o8", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clost5", "o__Clostridiales", "f__Lachnospiraceae"]}}],"columns": [{"id": "s3", "metadata": null},{"id": "s4", "metadata": null},{"id": "s5", "metadata": null}]}'
# BIOM_STRING_3 lacks samples 1 and 2 and has string metadata
BIOM_STRING_3 = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2-dev","date": "2013-04-17T20:23:50.481011","matrix_type": "sparse","matrix_element_type": "float","shape": [8, 3],"data": [[0,0,3.0],[0,1,4.0],[0,2,5.0],[1,0,8.0],[1,1,9.0],[1,2,10.0],[2,0,13.0],[2,1,14.0],[2,2,15.0],[3,0,18.0],[3,1,19.0],[3,2,20.0],[4,0,23.0],[4,1,24.0],[4,2,25.0],[5,0,28.0],[5,1,29.0],[5,2,30.0],[6,0,33.0],[6,1,34.0],[6,2,35.0],[7,0,38.0],[7,1,39.0],[7,2,40.0]],"rows": [{"id": "o1", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clost5;o__Clostridiales;f__Lachnospiraceae"}},{"id": "o2", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnos1"}},{"id": "o3", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnos1"}},{"id": "o4", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnos2"}},{"id": "o5", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnos2"}},{"id": "o6", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostri3;f__Lachnospiraceae"}},{"id": "o7", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostri3;f__Lachnospiraceae"}},{"id": "o8", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clost5;o__Clostridiales;f__Lachnospiraceae"}}],"columns": [{"id": "s3", "metadata": null},{"id": "s4", "metadata": null},{"id": "s5", "metadata": null}]}'
# BIOM_STRING_4 has metadata as dict
BIOM_STRING_4 = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "will","date": "2013-04-17T20:51:39.045533","matrix_type": "sparse","matrix_element_type": "float","shape": [8, 5],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,4.0],[0,4,5.0],[1,0,6.0],[1,1,7.0],[1,2,8.0],[1,3,9.0],[1,4,10.0],[2,0,11.0],[2,1,12.0],[2,2,13.0],[2,3,14.0],[2,4,15.0],[3,0,16.0],[3,1,17.0],[3,2,18.0],[3,3,19.0],[3,4,20.0],[4,0,21.0],[4,1,22.0],[4,2,23.0],[4,3,24.0],[4,4,25.0],[5,0,26.0],[5,1,27.0],[5,2,28.0],[5,3,29.0],[5,4,30.0],[6,0,31.0],[6,1,32.0],[6,2,33.0],[6,3,34.0],[6,4,35.0],[7,0,36.0],[7,1,37.0],[7,2,38.0],[7,3,39.0],[7,4,40.0]],"rows": [{"id": "o1", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clost5", "family": "f__Lachnospiraceae", "order": "o__Clostridiales"}}},{"id": "o2", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnos1", "order": "o__Clostridiales"}}},{"id": "o3", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnos1", "order": "o__Clostridiales"}}},{"id": "o4", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnos2", "order": "o__Clostridiales"}}},{"id": "o5", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnos2", "order": "o__Clostridiales"}}},{"id": "o6", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnospiraceae", "order": "o__Clostri3"}}},{"id": "o7", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnospiraceae", "order": "o__Clostri3"}}},{"id": "o8", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clost5", "family": "f__Lachnospiraceae", "order": "o__Clostridiales"}}}],"columns": [{"id": "s1", "metadata": null},{"id": "s2", "metadata": null},{"id": "s3", "metadata": null},{"id": "s4", "metadata": null},{"id": "s5", "metadata": null}]}'
# BIOM_STRING_5 has metadata as default dict though appears to be the same
# str??
BIOM_STRING_5 = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "will","date": "2013-04-17T21:20:27.650446","matrix_type": "sparse","matrix_element_type": "float","shape": [8, 5],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,4.0],[0,4,5.0],[1,0,6.0],[1,1,7.0],[1,2,8.0],[1,3,9.0],[1,4,10.0],[2,0,11.0],[2,1,12.0],[2,2,13.0],[2,3,14.0],[2,4,15.0],[3,0,16.0],[3,1,17.0],[3,2,18.0],[3,3,19.0],[3,4,20.0],[4,0,21.0],[4,1,22.0],[4,2,23.0],[4,3,24.0],[4,4,25.0],[5,0,26.0],[5,1,27.0],[5,2,28.0],[5,3,29.0],[5,4,30.0],[6,0,31.0],[6,1,32.0],[6,2,33.0],[6,3,34.0],[6,4,35.0],[7,0,36.0],[7,1,37.0],[7,2,38.0],[7,3,39.0],[7,4,40.0]],"rows": [{"id": "o1", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clost5", "family": "f__Lachnospiraceae", "order": "o__Clostridiales"}}},{"id": "o2", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnos1", "order": "o__Clostridiales"}}},{"id": "o3", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnos1", "order": "o__Clostridiales"}}},{"id": "o4", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnos2", "order": "o__Clostridiales"}}},{"id": "o5", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnos2", "order": "o__Clostridiales"}}},{"id": "o6", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnospiraceae", "order": "o__Clostri3"}}},{"id": "o7", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clostridia", "family": "f__Lachnospiraceae", "order": "o__Clostri3"}}},{"id": "o8", "metadata": {"taxonomy": {"kingdom": "k__Bacteria", "phylum": "p__Firmicutes", "class": "c__Clost5", "family": "f__Lachnospiraceae", "order": "o__Clostridiales"}}}],"columns": [{"id": "s1", "metadata": null},{"id": "s2", "metadata": null},{"id": "s3", "metadata": null},{"id": "s4", "metadata": null},{"id": "s5", "metadata": null}]}'
# BIOM_STRING_6 has s1 all zeros and o5 all zeros
BIOM_STRING_6 = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2-dev","date": "2013-04-17T23:22:41.269681","matrix_type": "sparse","matrix_element_type": "float","shape": [8, 5],"data": [[0,1,2.0],[0,2,3.0],[0,3,4.0],[0,4,5.0],[1,1,7.0],[1,2,8.0],[1,3,9.0],[1,4,10.0],[2,1,12.0],[2,2,13.0],[2,3,14.0],[2,4,15.0],[3,1,17.0],[3,2,18.0],[3,3,19.0],[3,4,20.0],[5,1,27.0],[5,2,28.0],[5,3,29.0],[5,4,30.0],[6,1,32.0],[6,2,33.0],[6,3,34.0],[6,4,35.0],[7,1,37.0],[7,2,38.0],[7,3,39.0],[7,4,40.0]],"rows": [{"id": "o1", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clost5;o__Clostridiales;f__Lachnospiraceae"}},{"id": "o2", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnos1"}},{"id": "o3", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnos1"}},{"id": "o4", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnos2"}},{"id": "o5", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnos2"}},{"id": "o6", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostri3;f__Lachnospiraceae"}},{"id": "o7", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostri3;f__Lachnospiraceae"}},{"id": "o8", "metadata": {"taxonomy": "k__Bacteria;p__Firmicutes;c__Clost5;o__Clostridiales;f__Lachnospiraceae"}}],"columns": [{"id": "s1", "metadata": null},{"id": "s2", "metadata": null},{"id": "s3", "metadata": null},{"id": "s4", "metadata": null},{"id": "s5", "metadata": null}]}'
MF_LINES = '#SampleID\tTreatment\tTimePt\tDiet\tStudy\ns1\tpre\t1\thf\ta\ns2\tpre\t2\tlf\ta\ns3\tpre\t3\thf\ta\ns4\tpost\t4\tlf\ta\ns5\tpost\t5\tmf\ta'


if __name__ == "__main__":
    main()
