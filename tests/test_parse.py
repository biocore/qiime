#!/usr/bin/env python
#file test_parse.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight", "Justin Kuczynski", "Greg Caporaso",\
                "Cathy Lozupone", "Jens Reeder"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"

from numpy import array, nan
from StringIO import StringIO
from cogent.util.unit_test import TestCase, main
from qiime.parse import (group_by_field, group_by_fields, 
    parse_distmat, parse_rarefaction_record, parse_rarefaction, parse_coords, 
    parse_otu_table, make_envs_dict, fields_to_dict, parse_rarefaction_fname,
    parse_qiime_parameters, parse_qiime_config_files,
    parse_bootstrap_support, parse_sample_mapping, parse_distmat_to_dict,
    sample_mapping_to_otu_table, parse_taxonomy, parse_mapping_file, 
    parse_metadata_state_descriptions, parse_rarefaction_data,
    parse_illumina_line, parse_qual_score, parse_qual_scores, QiimeParseError,
    parse_newick)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        self.l19_data = array([
            [7,1,0,0,0,0,0,0,0],
            [4,2,0,0,0,1,0,0,0],
            [2,4,0,0,0,1,0,0,0],
            [1,7,0,0,0,0,0,0,0],
            [0,8,0,0,0,0,0,0,0],
            [0,7,1,0,0,0,0,0,0],
            [0,4,2,0,0,0,2,0,0],
            [0,2,4,0,0,0,1,0,0],
            [0,1,7,0,0,0,0,0,0],
            [0,0,8,0,0,0,0,0,0],
            [0,0,7,1,0,0,0,0,0],
            [0,0,4,2,0,0,0,3,0],
            [0,0,2,4,0,0,0,1,0],
            [0,0,1,7,0,0,0,0,0],
            [0,0,0,8,0,0,0,0,0],
            [0,0,0,7,1,0,0,0,0],
            [0,0,0,4,2,0,0,0,4],
            [0,0,0,2,4,0,0,0,1],
            [0,0,0,1,7,0,0,0,0]
            ])
        self.l19_sample_names = ['sam1', 'sam2', 'sam3', 'sam4', 'sam5','sam6',\
        'sam7', 'sam8', 'sam9', 'sam_middle', 'sam11', 'sam12', 'sam13', \
        'sam14', 'sam15', 'sam16', 'sam17', 'sam18', 'sam19']
        self.l19_taxon_names =  ['tax1', 'tax2', 'tax3', 'tax4', 'endbigtaxon',\
        'tax6', 'tax7', 'tax8', 'tax9']
        self.SampleMapping = ["OTU1\tsample1\t3", "OTU1\tsample3\t2", \
        "OTU2\tsample1\t1", "OTU2\tsample2\t2"]
        self.SampleMapping2 = ["OTU1\tsample1", "OTU1\tsample3", \
        "OTU2\tsample1", "OTU2\tsample2"]
        
    def test_parse_newick(self):
        """parse_newick correctly matches escaped tip names to otu ids
        """
        # confirm that it works without escaped names
        t1 = ('((((tax7:0.1,tax3:0.2):.98,tax8:.3, tax4:.3):.4,'
              '((tax1:0.3, tax6:.09):0.43,tax2:0.4):0.5):.2,'
              '(tax9:0.3, endbigtaxon:.08));')
        expected1 = ['tax7','tax3','tax8','tax4','tax1',
                     'tax6','tax2','tax9','endbigtaxon']
        self.assertEqual(set(parse_newick(t1).getTipNames()),set(expected1))
        self.assertEqual(set([tip.Name for tip in parse_newick(t1).tips()]),
                         set(expected1))
                         
        # throw some screwed up names in 
        t2 = ('((((tax7:0.1,tax3:0.2):.98,tax8:.3, \'tax4\':.3):.4,'
              "(('ta_______ x1':0.3, tax6:.09):0.43,tax2:0.4):0.5):.2,"
              '(tax9:0.3, endbigtaxon:.08));')
        expected2 = ['tax7','tax3','tax8','tax4','ta_______ x1',
                     'tax6','tax2','tax9','endbigtaxon']
        self.assertEqual(set(parse_newick(t2).getTipNames()),set(expected2))
        self.assertEqual(set([tip.Name for tip in parse_newick(t2).tips()]),
                         set(expected2))
        

    def test_parse_mapping_file(self):
        """parse_mapping_file functions as expected"""
        s1 = ['#sample\ta\tb', '#comment line to skip',\
              'x \t y \t z ', ' ', '#more skip', 'i\tj\tk']
        exp = ([['x','y','z'],['i','j','k']],\
               ['sample','a','b'],\
               ['comment line to skip','more skip'])
        obs = parse_mapping_file(s1)
        self.assertEqual(obs, exp)
        
        #check that we strip double quotes by default
        s2 = ['#sample\ta\tb', '#comment line to skip',\
              '"x "\t" y "\t z ', ' ', '"#more skip"', 'i\t"j"\tk']
        obs = parse_mapping_file(s2)
        self.assertEqual(obs, exp)

    def test_group_by_field(self):
        """group_by_field should group table by fields"""
        t = [
                ['#sample', 'loc', 'age'],
                ['a','US','5'],
                ['b','US','10'],
                ['c','Mal','5'],
                ['d','Mal','10'],
                ['e','Ven','5'],
            ]
        self.assertEqual(group_by_field(t, 'loc'), \
            {'US':['a','b'], 'Mal':['c','d'], 'Ven':['e']})
        self.assertEqual(group_by_field(t, 'age'), \
            {'5':['a','c','e'], '10':['b','d']})

    def test_group_by_fields(self):
        """group_by_fields should group table by fields"""
        t = [
                ['#sample', 'loc', 'age', 'mal'],
                ['a','US','5','n'],
                ['b','US','10','n'],
                ['c','Mal','5','y'],
                ['d','Mal','10','n'],
                ['e','Mal','5','y'],
            ]
        self.assertEqual(group_by_fields(t, ['age','loc']), \
            {('5','US'):['a'], ('10','US'):['b'], ('5','Mal'):['c','e'],
            ('10','Mal'):['d']})

    def test_parse_distmat(self):
        """parse_distmat should read distmat correctly"""
        lines = """\ta\tb\tc
a\t0\t1\t2
b\t1\t0\t3.5
c\t1\t3.5\t0
""".splitlines()
        exp = (['a','b','c'], array([[0,1,2],[1,0,3.5],[1,3.5,0]]))
        obs = parse_distmat(lines)
        self.assertEqual(obs, exp)

    def test_parse_distmat_to_dict(self):
        """parse_distmat should return dict of distmat"""
        lines = """\ta\tb\tc
a\t0\t1\t2
b\t1\t0\t3.5
c\t1\t3.5\t0
""".splitlines()
        exp = {'a': {'a': 0.0, 'c': 2.0, 'b': 1.0},
                'c': {'a': 1.0, 'c': 0.0, 'b': 3.5},
                'b': {'a': 1.0, 'c': 3.5, 'b': 0.0}}
        obs = parse_distmat_to_dict(lines)
        self.assertEqual(obs, exp)

        #should raise error because row and column headers don't match
        wrong_dist_mat ="""\ta\ty\tx
a\t0\t1\t2
b\t1\t0\t3.5
c\t1\t3.5\t0
""".splitlines()
        self.failUnlessRaises(AssertionError, parse_distmat_to_dict, wrong_dist_mat)
        
    def test_parse_bootstrap_support(self):
        """parse_distmat should read distmat correctly"""
        input_txt = """#\ta\tb\tc.
#more comments here
node2\t0
17node\t0.11922
"""
        lines = input_txt.splitlines()
        exp = {'17node':0.11922, 'node2':0.00}
        obs = parse_bootstrap_support(lines)
        self.assertFloatEqual(obs, exp)

    def test_parse_rarefaction_data(self):
        self.data = {}
        self.data['headers'] = ['PD_whole_tree.txt', 'Antibiotics']
        self.data['error'] = {'NA': [0.099969643842700004], \
        'Y': [0.105669693476, 1.08546135424, 1.5626248357999999], \
        'N': [0.101173002974]}
        self.data['options'] = ['Y', 'NA', 'N']
        self.data['xaxis'] = [10.0, 310.0, 610.0, 910.0, 1210.0, 1510.0, \
        1810.0, 2110.0, 2410.0, 2710.0, 3010.0]
        self.data['series'] = {'NA': [0.88581050485400004], \
        'Y': [0.918845147059, 7.1758656176500004, 9.9186072941199992], \
        'N': [0.92636763785999998]}
        self.data['color'] = {'NA': '#00ff00', 'Y': '#ff0000', 'N': '#0000ff'}
        
        self.rarefaction_series_data = ['# PD_whole_tree.txt',
        '# Antibiotics',
        'xaxis: 10.0\t310.0\t610.0\t910.0\t1210.0\t1510.0\t1810.0\t2110.0\
        \t2410.0\t2710.0\t3010.0\t',
        'xmax: 3310.0',
        '>> Y',
        'color #ff0000',
        'series 0.918845147059\t7.17586561765\t9.91860729412\t',
        'error 0.105669693476\t1.08546135424\t1.5626248358\t',
        '>> NA',
        'color #00ff00',
        'series 0.885810504854\t',
        'error 0.0999696438427\t',
        '>> N',
        'color #0000ff',
        'series 0.92636763786\t',
        'error 0.101173002974'
        ]
        test = parse_rarefaction_data(self.rarefaction_series_data)
        self.assertEqual(test, self.data)

    def test_parse_rarefaction_record(self):
        self.rarefactionline1 = 'rare10.txt\t10\t0\t1.99181\t0.42877\t2.13996'
        test1 = parse_rarefaction_record(self.rarefactionline1)
        self.rarefactiondata1 = ('rare10.txt', [10.0, 0.0, \
        1.9918100000000001, 0.42876999999999998, 2.1399599999999999])
        self.assertEqual(self.rarefactiondata1, test1)
        
        self.rarefactionline2 = 'rare10.txt\t10\t0\t1.99181\t0.42877\tNA'
        test2 = parse_rarefaction_record(self.rarefactionline2)
        self.rarefactiondata2 = ('rare10.txt', [10.0, 0.0, 1.9918100000000001,\
         0.42876999999999998, nan])
        self.assertEqual(self.rarefactiondata2, test2)

    def test_parse_rarefaction_fname(self):
        """ parse_rarefaction_fname should return base, seqs/sam, iters, etc."""
        fname = "alpha_rarefaction_900_3.txt"
        base, seqs, iter, ext = parse_rarefaction_fname(fname)
        self.assertEqual((base, seqs, iter, ext),
            ("alpha_rarefaction", 900, 3, ".txt"))

    def test_parse_rarefaction(self):
        self.rarefactionfile = ['\tsequences per sample\titeration\t123\t234\t345',
                                'rare10.txt\t10\t0\t1.99181\t0.42877\t2.13996',
                                'rare10.txt\t10\t1\t2.07163\t0.42877\t2.37055',
                                'rare310.txt\t310\t0\t8.83115\t0.42877\t11.00725',
                                'rare310.txt\t310\t1\t10.05242\t0.42877\t8.24474',
                                'rare610.txt\t610\t0\t12.03067\t0.42877\t11.58928',
                                'rare610.txt\t610\t1\t12.9862\t0.42877\t11.58642']
        
        self.col_headers = ['', 'sequences per sample', 'iteration', '123', '234', '345']
        self.comments = []
        self.rarefaction_fns = ['rare10.txt', 'rare10.txt', 'rare310.txt', 'rare310.txt', 'rare610.txt', 'rare610.txt']
        self.rarefaction_data = [[10.0, 0.0, 1.9918100000000001, 0.42876999999999998, 2.1399599999999999], [10.0, 1.0, 2.0716299999999999, 0.42876999999999998, 2.3705500000000002], [310.0, 0.0, 8.8311499999999992, 0.42876999999999998, 11.007250000000001], [310.0, 1.0, 10.05242, 0.42876999999999998, 8.2447400000000002], [610.0, 0.0, 12.030670000000001, 0.42876999999999998, 11.58928], [610.0, 1.0, 12.9862, 0.42876999999999998, 11.58642]]
        
        test_col_headers, test_comments, test_rarefaction_fns, test_rarefaction_data = parse_rarefaction(self.rarefactionfile)
        self.assertEqual(test_col_headers, self.col_headers)
        self.assertEqual(test_comments, self.comments)
        self.assertEqual(test_rarefaction_fns, self.rarefaction_fns)
        self.assertEqual(test_rarefaction_data, self.rarefaction_data)

#     def test_parse_rarefaction(self):
#         """parse_rarefaction should handle multiple recs"""
#         recs ="""#HEADER  97.0    NFkeyRightShift 1000
# #CHAO1    288.12903   241.27093   371.41733
# #ACE  294.74813   252.75930   361.30606   0.68654
# #SHANNON  3.71822 3.61146 3.82499
# #SIMPSON  0.08021
# #n    rare    rare_lci    rare_hci
# 1 1.000000    1.000000    1.000000
# 51    26.878000   26.032290   27.723710
# #HEADER   97.0    DMkeySpace  1000
# #CHAO1    90.20000    53.81120    196.98011
# #ACE  122.66234   68.48901    264.46888   1.53571
# #SHANNON  1.35156 1.12549 1.57762
# #SIMPSON  0.56002
# #n    rare    rare_lci    rare_hci
# 1 1.000000    1.000000    1.000000
# 51    10.707000   10.085986   11.328014
# 101   17.547000   17.046632   18.047368
# 151   23.410000   23.009805   23.810195
# #HEADER   97.0    RKkeyW  1000
# #CHAO1    251.57143   176.86438   401.83780
# #ACE  264.44294   192.23364   395.08515   0.89714
# #SHANNON  1.56521 1.44411 1.68630
# #SIMPSON  0.54772
# #n    rare    rare_lci    rare_hci
# 1 1.000000    1.000000    1.000000
# 51    11.739000   11.060658   12.417342
# 101   19.054000   18.458137   19.649863
# 151   24.994000   24.447174   25.540826
# """.splitlines()
#         obs = parse_rarefaction(recs)
#         self.assertEqual(set(obs.keys()), \
#             set(['NFkeyRightShift','DMkeySpace','RKkeyW']))

    def test_parse_coords(self):
        """parse_coords should handle coords file"""
        coords = """pc vector number\t1\t2\t3
A\t0.11\t0.09\t0.23
B\t0.03\t0.07\t-0.26
C\t0.12\t0.06\t-0.32


eigvals\t4.94\t1.79\t1.50
% variation explained\t14.3\t5.2\t4.3


""".splitlines()
        obs = parse_coords(coords)
        exp = (['A','B','C'], 
            array([[.11,.09,.23],[.03,.07,-.26],[.12,.06,-.32]]),
            array([4.94,1.79,1.50]),
            array([14.3,5.2,4.3]))
        self.assertEqual(obs, exp)

    
    def test_parse_otu_table(self):
        """parse_otu_table should return correct result from small table"""
        data = """#Full OTU Counts
#OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors"""
        data_f = (data.split('\n'))
        obs = parse_otu_table(data_f)
        exp = (['Fing','Key','NA'],
               ['0','1','2','3','4'],
               array([[19111,44536,42],[1216,3500,6],[1803,1184,2],\
                    [1722,4903,17], [589,2074,34]]),
               [['Bacteria','Actinobacteria','Actinobacteridae','Propionibacterineae','Propionibacterium'],
                ['Bacteria','Firmicutes','Alicyclobacillaceae','Bacilli','Lactobacillales','Lactobacillales','Streptococcaceae','Streptococcus'],
                ['Bacteria','Actinobacteria','Actinobacteridae','Gordoniaceae','Corynebacteriaceae'],
                ['Bacteria','Firmicutes','Alicyclobacillaceae','Bacilli','Staphylococcaceae'],
                ['Bacteria','Cyanobacteria','Chloroplasts','vectors']])
        self.assertEqual(obs, exp)
        
    def test_parse_otu_table_float_counts(self):
        """parse_otu_table should return correct result from small table"""
        data = """#Full OTU Counts
#OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors"""
        data_f = (data.split('\n'))
        obs = parse_otu_table(data_f,count_map_f=float)
        exp = (['Fing','Key','NA'],
               ['0','1','2','3','4'],
               array([[19111.,44536.,42.],[1216.,3500.,6.],[1803.,1184.,2.],\
                    [1722.,4903.,17.], [589,2074.,34.]]),
               [['Bacteria','Actinobacteria','Actinobacteridae','Propionibacterineae','Propionibacterium'],
                ['Bacteria','Firmicutes','Alicyclobacillaceae','Bacilli','Lactobacillales','Lactobacillales','Streptococcaceae','Streptococcus'],
                ['Bacteria','Actinobacteria','Actinobacteridae','Gordoniaceae','Corynebacteriaceae'],
                ['Bacteria','Firmicutes','Alicyclobacillaceae','Bacilli','Staphylococcaceae'],
                ['Bacteria','Cyanobacteria','Chloroplasts','vectors']])
        self.assertEqual(obs, exp)
        
    def test_parse_otu_table_file(self):
        """parse_otu_table should return correct result on fileio format object"""
        data = """#Full OTU Counts
#OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors"""
        data_f = StringIO(data)
        obs = parse_otu_table(data_f)
        exp = (['Fing','Key','NA'],
               ['0','1','2','3','4'],
               array([[19111,44536,42],[1216,3500,6],[1803,1184,2],\
                    [1722,4903,17], [589,2074,34]]),
               [['Bacteria','Actinobacteria','Actinobacteridae','Propionibacterineae','Propionibacterium'],
                ['Bacteria','Firmicutes','Alicyclobacillaceae','Bacilli','Lactobacillales','Lactobacillales','Streptococcaceae','Streptococcus'],
                ['Bacteria','Actinobacteria','Actinobacteridae','Gordoniaceae','Corynebacteriaceae'],
                ['Bacteria','Firmicutes','Alicyclobacillaceae','Bacilli','Staphylococcaceae'],
                ['Bacteria','Cyanobacteria','Chloroplasts','vectors']])
        self.assertEqual(obs, exp)

    def test_make_envs_dict(self):
        """ make_envs_dict should have the same abundance for each taxon
        as the matrix that made the dict"""
        envs = make_envs_dict(self.l19_data, self.l19_sample_names,
            self.l19_taxon_names)
        for key in envs.keys():
            col_idx = self.l19_taxon_names.index(key)
            self.assertEqual(sum(envs[key].values()),
                self.l19_data[:,col_idx].sum())

    def test_fields_to_dict(self):
        """fields_to_dict should make first field key, rest val"""
        test_data = \
"""0	R27DLI_4812	R27DLI_600	R27DLI_727	U1PLI_403	U1PLI_8969	U1PLI_9080	U1PLI_9526	W3Cecum_6642	W3Cecum_8992
1	U1PLI_7889
2	W3Cecum_4858
3	R27DLI_3243	R27DLI_4562	R27DLI_6828	R27DLI_9097	U1PLI_2780	U1PLI_67	U9PSI_10475	U9PSI_4341	W3Cecum_5191""".splitlines()    #output from cd-hit
        obs = fields_to_dict(test_data)
        exp = {'0':['R27DLI_4812','R27DLI_600','R27DLI_727','U1PLI_403','U1PLI_8969','U1PLI_9080','U1PLI_9526','W3Cecum_6642','W3Cecum_8992'],
                '1':['U1PLI_7889'],
                '2':['W3Cecum_4858'],
                '3':['R27DLI_3243','R27DLI_4562','R27DLI_6828','R27DLI_9097','U1PLI_2780','U1PLI_67','U9PSI_10475','U9PSI_4341','W3Cecum_5191']}
        self.assertEqual(obs, exp)
        
    def test_parse_qiime_parameters(self):
        """parse_qiime_parameters: functions with valid input """
        lines = ["#Don't edit this file!",\
                 "pick_otus:similarity 0.94",\
                 "pick_otus:otu_picking_method\tcdhit",\
                 "align_seqs:verbose",\
                 "assign_taxonomy:use_rdp\ttRuE",\
                 "assign_taxonomy:something\tNone",\
                 "",\
                 "#some_script:fake_parameter\t99.0"]
        actual = parse_qiime_parameters(lines)
        expected = {'pick_otus':\
                     {'similarity':'0.94', 'otu_picking_method':'cdhit'},\
                    'assign_taxonomy':\
                     {'use_rdp':None}}
        self.assertEqual(actual,expected)
        
        # default dict functions as expected -- looking up non-existant key
        # returns empty dict
        self.assertEqual(actual['some_other_script'],{})

    def test_parse_sample_mapping(self):
        """parse_sample_mapping works"""
        lines = self.SampleMapping
        OTU_sample_info, all_sample_names = parse_sample_mapping(lines)
        self.assertEqual(OTU_sample_info, {'OTU2': {'sample1': '1', 'sample3': '0', 'sample2': '2'}, 'OTU1': {'sample1': '3', 'sample3': '2', 'sample2': '0'}})
        self.assertEqual(all_sample_names, set(['sample1', 'sample3', 'sample2']))
        #test that it works if no sample counts in there
        lines = self.SampleMapping2
        OTU_sample_info, all_sample_names = parse_sample_mapping(lines)
        self.assertEqual(OTU_sample_info, {'OTU2': {'sample1': '1', 'sample3': '0', 'sample2': '1'}, 'OTU1': {'sample1': '1', 'sample3': '1', 'sample2': '0'}})
        self.assertEqual(all_sample_names, set(['sample1', 'sample3', 'sample2']))
    def test_sample_mapping_to_otu_table(self):
        """sample_mapping_to_otu_table works"""
        lines = self.SampleMapping
        result = sample_mapping_to_otu_table(lines)
        self.assertEqual(result, ['#Full OTU Counts',\
         '#OTU ID\tsample1\tsample2\tsample3', 'OTU2\t1\t2\t0', \
        'OTU1\t3\t0\t2'])
    
    def test_parse_taxonomy(self):
        """ should parse taxonomy example, keeping otu id only"""
        example_tax = \
"""412 PC.635_647	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales	0.930
319 PC.355_281	Root;Bacteria;Bacteroidetes	0.970
353 PC.634_154	Root;Bacteria;Bacteroidetes	0.830
17 PC.607_302	Root;Bacteria;Bacteroidetes	0.960
13 PC.481_1214	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales	0.870
338 PC.593_1314	Root;Bacteria	0.990	42556	Additional fields ignored"""
        res = parse_taxonomy(example_tax.split('\n'))
        self.assertEqual(res['412'],
         "Root;Bacteria;Firmicutes;\"Clostridia\";Clostridiales")
        self.assertEqual(res['338'],
         "Root;Bacteria")
 
    def test_parse_qiime_config_files(self):
        """ parse_qiime_config_files functions as expected """
        fake_file1 = ['key1\tval1','key2 val2']
        fake_file2 = ['key2\tval3']
        actual = parse_qiime_config_files([fake_file1,fake_file2])
        expected = {'key1':'val1','key2':'val3'}
        self.assertEqual(actual,expected)
        
        # looking up a non-existant value returns None
        self.assertEqual(actual['fake_key'],None)
        
        # empty dict on empty input
        self.assertEqual(parse_qiime_config_files([]),{})


    def test_parse_metadata_state_descriptions(self):
        """parse_metadata_state_descriptions should return correct states from string."""
        s = ''
        self.assertEqual(parse_metadata_state_descriptions(s), {})
        s = 'Study:Twin,Hand,Dog;BodySite:Palm,Stool'
        self.assertEqual(parse_metadata_state_descriptions(s), {'Study':set(['Twin','Hand','Dog']),
            'BodySite':set(['Palm','Stool'])})
            
        
    def test_parse_illumina_line(self):
        """parse_illumina_line: functions with several lines """
        illumina_line0 = illumina_read1[0]
        illumina_line1 = illumina_read1[1]
        actual = parse_illumina_line(
         illumina_line0,barcode_length=6,rev_comp_barcode=True)
        expected = {\
         'Full description':'HWI-6X_9267:1:1:4:1699#ACCACCC/1',\
         'Machine Name':'HWI-6X_9267',\
         'Channel Number':1,\
         'Tile Number':1,\
         'X Position':4,\
         'Y Position':1699,\
         'Barcode':'GGTGGT',\
         'Full Y Position Field':'1699#ACCACCC/1',\
         'Sequence':\
          'TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCC'+\
          'CCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA',\
         'Quality Score':\
          'abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'+\
          'aaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB'}
        self.assertEqual(actual,expected)
        
        actual = parse_illumina_line(
         illumina_line0,barcode_length=6,rev_comp_barcode=False)
        expected['Barcode'] = 'ACCACC'
        
        actual = parse_illumina_line(
         illumina_line1,barcode_length=6,rev_comp_barcode=True)
        expected = {\
         'Full description':'HWI-6X_9267:1:1:4:390#ACCTCCC/1',\
         'Machine Name':'HWI-6X_9267',\
         'Channel Number':1,\
         'Tile Number':1,\
         'X Position':4,\
         'Y Position':390,\
         'Barcode':'GGAGGT',\
         'Full Y Position Field':'390#ACCTCCC/1',\
         'Sequence':\
          'GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGG'+\
          'GGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAA',\
         'Quality Score':\
          'aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'+\
          'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaBaaaaa'}
        self.assertEqual(actual,expected)
        
        actual = parse_illumina_line(
         illumina_line1,barcode_length=6,rev_comp_barcode=False)
        expected['Barcode'] = 'ACCTCC'

    def test_parse_qual_score(self):
        """qual_score should return dict of {id: qual_scores}"""
        scores = StringIO('>x\n5 10 5\n12\n>y\n30 40')
        self.assertEqual(parse_qual_score(scores),
                         {'x':[5,10,5,12],'y':[30,40]})

        #Check that a bad file, e.g. a fast raises Error
        bad_scores = StringIO('>x\nabcbd\n12\n>y\GATC')
        self.assertRaises(QiimeParseError, parse_qual_score, bad_scores)

    def test_parse_qual_scores(self):
        """qual_scores should return dict of {id:qual_scores}"""
        scores = StringIO('>x\n5 10 5\n12\n>y\n30 40')
        scores2= StringIO('>a\n5 10 5\n12\n>b\n30 40')
        self.assertEqual(parse_qual_scores([scores, scores2]),
            {'x':[5,10,5,12],'y':[30,40],'a':[5,10,5,12],'b':[30,40]})


illumina_read1 = """HWI-6X_9267:1:1:4:1699#ACCACCC/1:TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA:abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB
HWI-6X_9267:1:1:4:390#ACCTCCC/1:GACAGGAGGAGCAAGTGTTATTCAAATTATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAA:aaaaaaaaaa```aa\^_aa``aVaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaBaaaaa""".split('\n')

illumina_read2 = """HWI-6X_9267:1:1:4:1699#ACCACCC/2:TTTTAAAAAAAAGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTAAAAAAAAACCCCCCCGGGGGGGGTTTTTTTAATTATTC:aaaaaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbcccccccccccccccccBcccccccccccccccc```````BBBB
HWI-6X_9267:1:1:4:390#ACCTCCC/2:ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG:aaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbb""".split('\n')


if __name__ =='__main__':
    main()
