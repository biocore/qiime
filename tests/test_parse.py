#!/usr/bin/env python
#file test_parse.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight", "Justin Kuczynski", "Greg Caporaso",\
                "Cathy Lozupone", "Jens Reeder"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.91"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"

from numpy import array
from StringIO import StringIO
from cogent.util.unit_test import TestCase, main
from qiime.parse import (parse_map, group_by_field, group_by_fields, 
    parse_distmat, parse_rarefaction_rec, parse_rarefaction, parse_coords, 
    otu_file_to_lineages, parse_otus, otu_table_to_envs, parse_sequences_by_otu,
    make_envs_dict, fields_to_dict, parse_rarefaction_fname, envs_to_otu_counts,
    otu_counts_to_matrix, envs_to_matrix, parse_qiime_parameters, 
    parse_bootstrap_support, parse_sample_mapping, parse_distmat_to_dict,
    sample_mapping_to_otu_table, parse_taxonomy)

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

    def test_parse_map(self):
        """parse_map should handle tab-delimited file as expected."""
        s = """#sample\ta\tb
#comment line to skip
x \t y \t z 

#more skip
i\tj\tk"""
        s2= """#sample\ta\tb
#comment line to skip
"x "\t" y "\t z 

"#more skip"
i\t"j"\tk"""
        exp = [['#sample','a','b'],['x','y','z'],['i','j','k']]
        obs = parse_map(s.splitlines())
        self.assertEqual(obs, exp)
        exp2 = ([['#sample','a','b'],['x','y','z'],['i','j','k']], \
                ['comment line to skip','more skip'])
        obs = parse_map(s.splitlines(), return_header=True)
        self.assertEqual(obs, exp2)
        #check that we strip double quotes by default
        obs = parse_map(s2.splitlines(), return_header=True)
        self.assertEqual(obs, exp2)

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

    def test_parse_rarefaction_rec(self):
        """parse_rarefaction_rec should produce expected results"""
        rec = """#HEADER	97.0	NFkeyRightShift	1000
#CHAO1	288.12903	241.27093	371.41733
#ACE	294.74813	252.75930	361.30606	0.68654
#SHANNON	3.71822	3.61146	3.82499
#SIMPSON	0.08021
#n	rare	rare_lci	rare_hci
    1	1.000000	1.000000	1.000000
    51	26.878000	26.032290	27.723710
    101	44.007000	43.212590	44.801410""".splitlines()
        exp = {'pct_sim':97.0,'sample_id':'NFkeyRightShift','num_iters':1000,
            'ACE':[294.74813,252.75930,361.30606,0.68654],
            'CHAO1':[288.12903,241.27093,371.41733],
            'SHANNON':[3.71822,3.61146,3.82499],
            'SIMPSON':[0.08021],
            'rarefaction_data':array([[1,1,1,1],[51,26.878000,26.032290,27.723710],
                [101,44.007000,43.212590,44.801410]])}
        obs = parse_rarefaction_rec(rec)
        self.assertEqual(set(obs.keys()), set(exp.keys()))
        for k, v in obs.items():
            self.assertFloatEqual(v, exp[k])

    def test_parse_rarefaction_fname(self):
        """ parse_rarefaction_fname should return base, seqs/sam, iters, etc."""
        fname = "alpha_rarefaction_900_3.txt"
        base, seqs, iter, ext = parse_rarefaction_fname(fname)
        self.assertEqual((base, seqs, iter, ext),
            ("alpha_rarefaction", 900, 3, ".txt"))

    def test_parse_rarefaction(self):
        """parse_rarefaction should handle multiple recs"""
        recs ="""#HEADER	97.0	NFkeyRightShift	1000
#CHAO1	288.12903	241.27093	371.41733
#ACE	294.74813	252.75930	361.30606	0.68654
#SHANNON	3.71822	3.61146	3.82499
#SIMPSON	0.08021
#n	rare	rare_lci	rare_hci
1	1.000000	1.000000	1.000000
51	26.878000	26.032290	27.723710
#HEADER	97.0	DMkeySpace	1000
#CHAO1	90.20000	53.81120	196.98011
#ACE	122.66234	68.48901	264.46888	1.53571
#SHANNON	1.35156	1.12549	1.57762
#SIMPSON	0.56002
#n	rare	rare_lci	rare_hci
1	1.000000	1.000000	1.000000
51	10.707000	10.085986	11.328014
101	17.547000	17.046632	18.047368
151	23.410000	23.009805	23.810195
#HEADER	97.0	RKkeyW	1000
#CHAO1	251.57143	176.86438	401.83780
#ACE	264.44294	192.23364	395.08515	0.89714
#SHANNON	1.56521	1.44411	1.68630
#SIMPSON	0.54772
#n	rare	rare_lci	rare_hci
1	1.000000	1.000000	1.000000
51	11.739000	11.060658	12.417342
101	19.054000	18.458137	19.649863
151	24.994000	24.447174	25.540826
""".splitlines()
        obs = parse_rarefaction(recs)
        self.assertEqual(set(obs.keys()), \
            set(['NFkeyRightShift','DMkeySpace','RKkeyW']))

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

    def test_otu_file_to_lineages(self):
        """otu_file_to_lineages should extract correct info"""
        data = """##################
# OTU 0 section start
# 100.00% from Bacteria; Proteobacteria; Gammaproteobacteria 
# 1 total lineages 
# 1 unique seqs in OTU
# -- Lineage Detail ---
# 1 orig seqs in OTU
# -- IDs for Arb ---
# (6140)
# 1 (100.00%) Bacteria; Proteobacteria; Gammaproteobacteria
##################
>AK197_6140 GI:-1|EVAL:-1|BITS:-1.0|ID:-1.0%|COV:-1.0%|LINS:Bacteria; Proteobacteria; Gammaproteobacteria
CTGGGCCGTGTCTCAGTCCCAGTGTGGCTGATCGTCCTCTCAGACCAGCTACGGATCGTCGCCTTGGTAGGCCGTTACCCCACCAACAAGCTAATCCGCCGCGGGCTCATCCTCGGTCGGAGGCCGAAGCTACCTTTCCCTTCGAGACCCGAAGGTCCAAAGGGCCATTCCGTATTAATCCGGGTTTCCCCGGGCTATCCGGATACCAAGGGCAGATTACCCACGTGTTACTCACCCGTTCGCCGCTTTCCCCGGTCCCGAAGGG
##################
# OTU -1 section start
# 100.00% from Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales 
# 2 total lineages 
# 3 unique seqs in OTU
# -- Lineage Detail ---
# 7 orig seqs in OTU
# -- IDs for Arb ---
# (574, 9547, 11347, 6452, 16811, 12887, 7067)
# 6 (85.71%) Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales
# 1 (14.29%) Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Actinomycineae; Actinomycetaceae
##################
>AK595_574 GI:-1|EVAL:-1|BITS:-1.0|ID:-1.0%|COV:-1.0%|LINS:Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales
CTGGGCCGTATCTCAGTCCCAGTGTGACCGAACACCCTCTCAGGCCGGTTACCCGTCCACGCCTTGGTAGGCCATCACCCCACCAACAAACTGATAGGCCGCAAGCCCATCCCCCACCAACCCCCAAAAGGAAGCCTTTCCAAACCCCACCATGCAGCAGGAAAAGAATATTCGGTATTAGCCCACGTTTCCGCGAGTTATCCCAAAGATGAAGGGCAGGTTACTACGTGTTACTCACCGTTCGCACTAGAA
>AK894_9547 GI:-1|EVAL:-1|BITS:-1.0|ID:-1.0%|COV:-1.0%|LINS:Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Actinomycineae; Actinomycetaceae
CTGGGCCGTATCTCAGTCCCAGTGTGACCGAACACCCTCTCAGGCCGGTTACCCGTCCACGCCTTGGTAGGCCATCACCCCACCAACAAACTGATAGGCCGCAAGCCCATCCCCCACCAACCCCCAAAAAGGAGGCCTTTCCAAACCCCACCATGCAGCAGGAAAAGAATATCCCGTATTAGCCCACGTTTCCGCGAGTTATCCAGAAGATGAAGGGCAGGTTACTTACGTGTTACTCACCCGTTCGCCAA
>AK7929_11347@@AK1194_6452@@AK33210_16811@@AK1194_12887@@AK3410_7067 GI:-1|EVAL:-1|BITS:-1.0|ID:-1.0%|COV:-1.0%|LINS:Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales
CTGGGCCGTATCTCAGTCCCAGTGTGACCGAACACCCTCTCAGGCCGGTTACCCGTCCACGCCTTGGTAGGCCATCACCCCACCAACAAACTGATAGGCCGCAAGCCCATCCCCCACCAACCCCCAAAAGGAAGCCTTTCCAAACCCCACCATGCAGCAGGAAAAGAATATTCGGTATTAGCCCACGTTTCCGCGAGTTATCCCAAAGATGAAGGGCAGGTTACTTACGTGTTACTCACCCGTTCGCCACTAGAAACGCCCCCGG
##################
# OTU -6150 section start
# 100.00% from Unclassified-Screened 
# 1 total lineages 
# 1 unique seqs in OTU
# -- Lineage Detail ---
# 1 orig seqs in OTU
# -- IDs for Arb ---
# (7266)
# 1 (100.00%) Unclassified-Screened
##################
>AK1194_7266 GI:-1|EVAL:-1|BITS:-1.0|ID:-1.0%|COV:-1.0%|LINS:Unclassified-Screened
CTGGCCCGTGTTTCAGTGCCAGTGTGGCCGGTCGCCCTCTCAGGCCGGCTACTGATCGATGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATCAGCC"""
        data_f = StringIO(data)
        obs = otu_file_to_lineages(data_f)
        exp = {'0':[['Bacteria','Proteobacteria','Gammaproteobacteria'],100.0],
               '-1': [['Bacteria', 'Actinobacteria', 'Actinobacteria', \
                       'Actinobacteridae','Actinomycetales'], 100.0],
               '-6150': [['Unclassified-Screened'], 100.0]}
        self.assertEqual(obs.keys(), exp.keys())
        self.assertEqual(obs, exp)
    
    def test_parse_otus(self):
        """parse_otus should return correct result from small table"""
        data = """#Full OTU Counts
#OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors"""
        data_f = (data.split('\n'))
        obs = parse_otus(data_f)
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
        
    def test_parse_otus_file(self):
        """parse_otus should return correct result on fileio format object"""
        data = """#Full OTU Counts
#OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
2	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors"""
        data_f = StringIO(data)
        obs = parse_otus(data_f)
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

    def test_otu_table_to_envs(self):
        """otu_table_to_envs should produce correct dict for fast_unifrac"""
        obs = otu_table_to_envs(['a','b','c'],['1','2'], \
            array([[0,2,3],[4,0,0]]))
        exp = {'1':{'b':2,'c':3},'2':{'a':4}}
        self.assertEqual(obs, exp)

    def test_parse_sequences_by_otu(self):
        """parse_sequences_by_otu should return expected results"""
        data = """##################
# OTU -29479 section start
# 100.00% from Bacteria; Firmicutes; Clostridia; Clostridiales; Ruminococcaceae; Anaerofilum 
# 1 total lineages 
# 1 unique seqs in OTU
# -- Lineage Detail ---
# 2 orig seqs in OTU
# -- IDs for Arb ---
# (489870, 510240)
# 2 (100.00%) Bacteria; Firmicutes; Clostridia; Clostridiales; Ruminococcaceae; Anaerofilum
##################
>M21Fotl_489870@@M21Fotl_510240 GI:-1|EVAL:-1|BITS:-1.0|ID:-1.0%|COV:-1.0%|LINS:Bacteria; Firmicutes; Clostridia; Clostridiales; Ruminococcaceae; Anaerofilum
CTGGGCCGTATCTCAGTCCCAATGTGGCCGATCAACCTCTCAGTCCGGCTACTGATCGTCGCCATGGTGGGCCGTTACCCCGCCATCTAGCTAATCAGACGCGAGCCCATCTCAGAGCACATAAAGCTTTGATCTTCAGAAGATGCCTTCCGAAGATGTTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCACT
##################
# OTU -29469 section start
# 100.00% from Bacteria; Fusobacteria; Fusobacteria; Fusobacteriales; Fusobacteriaceae; Fusobacterium 
# 1 total lineages 
# 3 unique seqs in OTU
# -- Lineage Detail ---
# 4 orig seqs in OTU
# -- IDs for Arb ---
# (403492, 60020, 237665, 270713)
# 4 (100.00%) Bacteria; Fusobacteria; Fusobacteria; Fusobacteriales; Fusobacteriaceae; Fusobacterium
##################
>F32Indr_403492 GI:-1|EVAL:-1|BITS:-1.0|ID:-1.0%|COV:-1.0%|LINS:Bacteria; Fusobacteria; Fusobacteria; Fusobacteriales; Fusobacteriaceae; Fusobacterium
CTGGGCCGTGTCTCAGTCCCAGTGTGGCTGATCACCCTCTCAGGCCGGCTACCCATCATCGCCTTGGTGAGCCGTTACCTCTCCAACTAGCTAATGGGACGCAAAGCTCTCTCACAGTGCATATAGCTTTCATAATCTTAGGATGCCCTAAAATCATAATATCAGGTATTAGCATTCGTTTCCAAATGTTGTCCCTGTCT
>M54Plmr_60020@@M54Plmr_237665 GI:-1|EVAL:-1|BITS:-1.0|ID:-1.0%|COV:-1.0%|LINS:Bacteria; Fusobacteria; Fusobacteria; Fusobacteriales; Fusobacteriaceae; Fusobacterium
CTGGACCGTGTCTCAGTTCCAGTGTGGCCGATCACCCTCTCAGGCCGGCTACCCATCATCGCCTTGGTGAGCCGTTACCTCTCCAACTAGCTAATGGGACGCAAAGCTCTCTCACAGCGCATATAGCTTTCATAATCCTAGGATGCCCTAAAATCATAATATCAGGTATTAGCATTCGTTTCCAAATGTTGTCCCTATCT
>F32Mout_270713 GI:-1|EVAL:-1|BITS:-1.0|ID:-1.0%|COV:-1.0%|LINS:Bacteria; Fusobacteria; Fusobacteria; Fusobacteriales; Fusobacteriaceae; Fusobacterium
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACCCATCATCGCCTTGGTGAGCCGTTACCTCTCCAACTAGCTAATGGGACGCAAAGCTCTCTCACAGCGCATATAGCTTTCATAATCTTAGGATGCCCTAAAATCATAATATCAGGTATTAGCATTCGTTTCCAAATGTTGTCCCTATCT
"""
        result = parse_sequences_by_otu(data.splitlines())
        self.assertEqual(result[0], {-29479:['M21Fotl_489870','M21Fotl_510240'],
           -29469:['F32Indr_403492','M54Plmr_60020','M54Plmr_237665','F32Mout_270713']})
        self.assertEqual(result[1], {'M21Fotl_489870':-29479,'M21Fotl_510240':-29479,\
            'F32Indr_403492':-29469,'M54Plmr_60020':-29469,'M54Plmr_237665':-29469,\
            'F32Mout_270713':-29469})


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


    def test_envs_to_otu_counts(self):
        """envs_to_otu_counts should produce right dict"""
        s="""s01\ta\t3
s02\ta\t1
s01\tb\t4
s03\tc\t5""".splitlines()
        res = envs_to_otu_counts(s)
        self.assertEqual(res, {('a','s01'):3,('a','s02'):1,('b','s01'):4,
            ('c','s03'):5})

    def test_otu_counts_to_matrix(self):
        """otu_counts_to_matrix should produce right matrix/order"""
        data = {('a','s01'):3,('a','s02'):1,('b','s01'):4, ('c','s03'):5}
        matrix, all_otus, all_sampleids = otu_counts_to_matrix(data)
        self.assertEqual(all_sampleids, ['a','b','c'])
        self.assertEqual(all_otus, ['s01','s02','s03'])
        self.assertEqual(matrix, array([[3,4,0],[1,0,0],[0,0,5]]))

    def test_envs_to_matrix(self):
        """envs_to_matrix should take envs file, convert to OTU matrix"""
        s="""s01\ta\t3
s02\ta\t1
s01\tb\t4
s03\tc\t5""".splitlines()
        matrix, all_otus, all_sampleids = envs_to_matrix(s)
        self.assertEqual(all_sampleids, ['a','b','c'])
        self.assertEqual(all_otus, ['s01','s02','s03'])
        self.assertEqual(matrix, array([[3,4,0],[1,0,0],[0,0,5]]))
        
    def test_parse_qiime_parameters(self):
        """parse_qiime_parameters: functions with valid input """
        lines = ["#Don't edit this file!",\
                 "pick_otus:similarity\t0.94",\
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

    def test_parse_sample_mapping(self):
        """parse_sample_mapping works"""
        lines = self.SampleMapping
        OTU_sample_info, all_sample_names = parse_sample_mapping(lines)
        self.assertEqual(OTU_sample_info, {'OTU2': {'sample1': '1', 'sample3': '0', 'sample2': '2'}, 'OTU1': {'sample1': '3', 'sample3': '2', 'sample2': '0'}})

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
338 PC.593_1314	Root;Bacteria	0.990"""
        res = parse_taxonomy(example_tax.split('\n'))
        self.assertEqual(res['412'],
         "Root;Bacteria;Firmicutes;\"Clostridia\";Clostridiales")
        self.assertEqual(res['338'],
         "Root;Bacteria")
 
if __name__ =='__main__':
    main()
