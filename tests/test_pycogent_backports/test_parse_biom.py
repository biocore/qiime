#!/usr/bin/env python
#file test_parse.py

__author__ = "Justin Kuczynski"
__copyright__ = ""
__credits__ = ["Justin Kuczynski"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

from numpy import array, nan
from StringIO import StringIO
from cogent.util.unit_test import TestCase,main
from cogent.util.misc import remove_files
from qiime.util import get_tmp_filename
from qiime.parse import parse_otu_table
import qiime.pycogent_backports.parse_biom as parse
import qiime.pycogent_backports.rich_otu_table as rich_otu_table

class ParseTests(TestCase):
    """Tests of parse functions"""

    def setUp(self):
        """define some top-level data"""
        self.legacy_otu_table1 = legacy_otu_table1
        self.otu_table1 = otu_table1
        self.otu_table1_floats=otu_table1_floats
        self.files_to_remove = []
        self.biom_minimal_dense = biom_minimal_dense
        self.biom_minimal_sparse = biom_minimal_sparse
        
        self.classic_otu_table1_w_tax = classic_otu_table1_w_tax.split('\n')
        self.classic_otu_table1_no_tax = classic_otu_table1_no_tax.split('\n')
        self.biom_otu_table1_w_tax = biom_otu_table1_w_tax.split('\n')
        self.biom_otu_table1_no_tax = biom_otu_table1_no_tax.split('\n')
    
    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_parse_otu_table_to_rich_otu_table(self):
        tab1_fh = StringIO(self.otu_table1)
        sparse_rich = parse.parse_otu_table_to_rich_otu_table(tab1_fh)
        self.assertEqual(sorted(sparse_rich.SampleIds),sorted(['Fing','Key','NA']))
        self.assertEqual(sorted(sparse_rich.ObservationIds),map(str,[0,1,3,4,7]))
        for i, obs_id in enumerate(sparse_rich.ObservationIds):
            if obs_id == '0':
                self.assertEqual(sparse_rich.ObservationMetadata[i],
                {'taxonomy':'Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium'.split('; ')})
            elif obs_id == '1':
                self.assertEqual(sparse_rich.ObservationMetadata[i],
                {'taxonomy':'Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus'.split('; ')})
            elif obs_id == '7':
                self.assertEqual(sparse_rich.ObservationMetadata[i],
                {'taxonomy':'Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae'.split('; ')})
            elif obs_id in ['3','4']: 
                pass # got lazy
            else: 
                raise RuntimeError('obs_id incorrect?')

        self.assertEquals(sparse_rich.SampleMetadata,None)

        for i, obs_id in enumerate(sparse_rich.ObservationIds):
            for j, sample_id in enumerate(sparse_rich.SampleIds):
                if obs_id == '1' and sample_id == 'Key':
                    self.assertEqual(True,True) # should test some abundance data

    def test_parse_otu_table_to_rich_otu_table_dense(self):
        tab1_fh = StringIO(self.otu_table1)
        sparse_rich = parse.parse_otu_table_to_rich_otu_table(tab1_fh,dense=True)
        self.assertEqual(sorted(sparse_rich.SampleIds),sorted(['Fing','Key','NA']))
        self.assertEqual(sorted(sparse_rich.ObservationIds),map(str,[0,1,3,4,7]))
        for i, obs_id in enumerate(sparse_rich.ObservationIds):
            if obs_id == '0':
                self.assertEqual(sparse_rich.ObservationMetadata[i],
                {'taxonomy':'Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium'.split('; ')})
            elif obs_id == '1':
                self.assertEqual(sparse_rich.ObservationMetadata[i],
                {'taxonomy':'Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus'.split('; ')})
            elif obs_id == '7':
                self.assertEqual(sparse_rich.ObservationMetadata[i],
                {'taxonomy':'Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae'.split('; ')})
            elif obs_id in ['3','4']: 
                pass # got lazy
            else: 
                raise RuntimeError('obs_id incorrect?')

        self.assertEquals(sparse_rich.SampleMetadata,None)

        for i, obs_id in enumerate(sparse_rich.ObservationIds):
            for j, sample_id in enumerate(sparse_rich.SampleIds):
                if obs_id == '1' and sample_id == 'Key':
                    self.assertEqual(True,True) # should test some abundance data


    def test_parse_biom_table_minimal_dense(self):
        tab1_fh = StringIO(self.biom_minimal_dense)
        tab = parse.parse_biom_table(tab1_fh)
        self.assertEqual((tab.SampleIds),['Sample1','Sample2',
            'Sample3','Sample4','Sample5','Sample6',])
        self.assertEqual((tab.ObservationIds),['GG_OTU_1','GG_OTU_2',
            'GG_OTU_3','GG_OTU_4','GG_OTU_5'])
        self.assertEqual(tab.SampleMetadata,None)
        self.assertEqual(tab.ObservationMetadata,None)

    def test_parse_biom_table_compare_sparse_dense(self):
        tab1_fh = StringIO(self.biom_minimal_dense)
        tab1 = parse.parse_biom_table(tab1_fh)
        tab2_fh = StringIO(self.biom_minimal_sparse)
        tab2 = parse.parse_biom_table(tab2_fh)
        assert(isinstance(tab1,rich_otu_table.DenseOTUTable))
        assert(isinstance(tab2,rich_otu_table.SparseOTUTable))

        tab1_sams = [sam[0] for sam in tab1.iterSamples()]
        tab2_sams = [sam[0] for sam in tab2.iterSamples()]

        self.assertFloatEqual(tab1_sams,tab2_sams)

        self.assertEqual((tab1.SampleIds),['Sample1','Sample2',
            'Sample3','Sample4','Sample5','Sample6',])
        self.assertEqual((tab1.ObservationIds),['GG_OTU_1','GG_OTU_2',
            'GG_OTU_3','GG_OTU_4','GG_OTU_5'])
        self.assertEqual((tab2.SampleIds),['Sample1','Sample2',
            'Sample3','Sample4','Sample5','Sample6',])
        self.assertEqual((tab2.ObservationIds),['GG_OTU_1','GG_OTU_2',
            'GG_OTU_3','GG_OTU_4','GG_OTU_5'])

    def test_otu_table_biom_conversions_with_taxonomy(self):
        """ converting between classic otu table and biom is roundtrip-able (w taxonomy)
        """
        # parse the classic otu table (w tax)
        parsed_classic_otu_table1_w_tax = parse_otu_table(self.classic_otu_table1_w_tax)
        
        # convert the classic otu table to a biom file, and then convert the biom file
        # to a classic otu table (i.e., roundtrip the file)
        roundtripped_classic_otu_table1_w_tax =\
          parse.convert_biom_to_otu_table(
            parse.convert_otu_table_to_biom(self.classic_otu_table1_w_tax))
        
        # parse the roundtripped file
        parsed_roundtripped_classic_otu_table1_w_tax = \
         parse_otu_table(roundtripped_classic_otu_table1_w_tax.split('\n'))
        
        # compare the parsed roundtripped file to the parsed input file
        self.assertEqual(parsed_roundtripped_classic_otu_table1_w_tax,
                         parsed_classic_otu_table1_w_tax)

    def test_otu_table_biom_conversions_no_taxonomy(self):
        """ converting between classic otu table and biom is roundtrip-able (no taxonomy)
        """
        # parse the classic otu table (no tax)
        parsed_classic_otu_table1_no_tax = parse_otu_table(self.classic_otu_table1_no_tax)
        
        # convert the classic otu table to a biom file, and then convert the biom file
        # to a classic otu table (i.e., roundtrip the file)
        roundtripped_classic_otu_table1_no_tax =\
          parse.convert_biom_to_otu_table(
            parse.convert_otu_table_to_biom(self.classic_otu_table1_no_tax))
        # parse the roundtripped file
        parsed_roundtripped_classic_otu_table1_no_tax = \
         parse_otu_table(roundtripped_classic_otu_table1_no_tax.split('\n'))
        
        # compare the parsed roundtripped file to the parsed input file
        self.assertEqual(parsed_roundtripped_classic_otu_table1_no_tax,
                         parsed_classic_otu_table1_no_tax)

legacy_otu_table1 = """# some comment goes here
#OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium

1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
7	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

otu_table1 = """# Some comment




OTU ID	Fing	Key	NA	Consensus Lineage
0	19111	44536	42	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
# some other comment
1	1216	3500	6	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
7	1803	1184	2	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
# comments
#    everywhere!
3	1722	4903	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589	2074	34	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

otu_table1_floats = """# Some comment




OTU ID	Fing	Key	NA	Consensus Lineage
0	19111.0	44536.0	42.0	Bacteria; Actinobacteria; Actinobacteridae; Propionibacterineae; Propionibacterium
# some other comment
1	1216.0	3500.0	6.0	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Lactobacillales; Lactobacillales; Streptococcaceae; Streptococcus
7	1803.0	1184.0	2.0	Bacteria; Actinobacteria; Actinobacteridae; Gordoniaceae; Corynebacteriaceae
# comments
#    everywhere!
3	1722.1	4903.2	17	Bacteria; Firmicutes; Alicyclobacillaceae; Bacilli; Staphylococcaceae
4	589.6	2074.4	34.5	Bacteria; Cyanobacteria; Chloroplasts; vectors
"""

biom_minimal_dense = """

    {
        "id":null,
        "format": "Biological Observation Matrix v0.9",
        "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "type": "OTU table",
        "generated_by": "QIIME revision XYZ",
        "date": "2011-12-19T19:00:00",
        "rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
            ],  
        "columns": [
                {"id":"Sample1", "metadata":null},
                {"id":"Sample2", "metadata":null},
                {"id":"Sample3", "metadata":null},
                {"id":"Sample4", "metadata":null},
                {"id":"Sample5", "metadata":null},
                {"id":"Sample6", "metadata":null}
            ],  
        "matrix_type": "dense",
        "matrix_element_type": "int",
        "shape": [5,6],
        "data":  [[0,0,1,0,0,0], 
                  [5,1,0,2,3,1],
                  [0,0,1,4,2,0],
                  [2,1,1,0,0,1],
                  [0,1,1,0,0,0]]
    }
"""

biom_minimal_sparse="""
    {
        "id":null,
        "format": "Biological Observation Matrix v0.9",
        "format_url": "http://some_website/QIIME_MGRAST_dataformat_v0.9.html",
        "type": "OTU table",
        "generated_by": "QIIME revision XYZ",
        "date": "2011-12-19T19:00:00",
        "rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
            ],  
        "columns": [
                {"id":"Sample1", "metadata":null},
                {"id":"Sample2", "metadata":null},
                {"id":"Sample3", "metadata":null},
                {"id":"Sample4", "metadata":null},
                {"id":"Sample5", "metadata":null},
                {"id":"Sample6", "metadata":null}
            ],
        "matrix_type": "sparse",
        "matrix_element_type": "int",
        "shape": [5, 6], 
        "data":[[0,2,1],
                [1,0,5],
                [1,1,1],
                [1,3,2],
                [1,4,3],
                [1,5,1],
                [2,2,1],
                [2,3,4],
                [2,4,2],
                [3,0,2],
                [3,1,1],
                [3,2,1],
                [3,5,1],
                [4,1,1],
                [4,2,1]
               ]
    }
"""

classic_otu_table1_w_tax = """#Full OTU Counts
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
0	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
1	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
2	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
3	2	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
4	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
5	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
6	0	0	0	0	0	0	0	1	0	Root;Bacteria;Actinobacteria;Actinobacteria
7	0	0	2	0	0	0	0	0	2	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
8	1	1	0	2	4	0	0	0	0	Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus
9	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
10	0	1	0	0	0	0	0	0	0	Root;Bacteria
11	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
12	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
13	1	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
14	0	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
15	0	0	0	0	1	0	0	0	0	Root;Bacteria
16	1	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
17	0	0	0	1	0	0	4	10	37	Root;Bacteria;Bacteroidetes
18	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
19	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
20	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
21	0	0	0	0	0	0	2	3	2	Root;Bacteria;Bacteroidetes
22	0	0	0	0	2	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
23	14	1	14	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus
24	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
25	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
26	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
27	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
28	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
29	6	0	4	0	2	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
30	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
31	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
32	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
33	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
34	0	0	0	0	0	0	8	10	2	Root;Bacteria
35	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
36	1	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
37	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
38	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
39	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
40	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
41	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
42	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
43	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
44	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
45	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales;Erysipelotrichaceae;Coprobacillus
46	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
47	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
48	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
49	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
50	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
51	0	1	0	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
52	0	2	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
53	0	0	0	0	0	0	2	0	1	Root;Bacteria;Proteobacteria;Deltaproteobacteria
54	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
55	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
56	0	0	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes
57	0	0	0	0	0	0	0	1	0	Root;Bacteria
58	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
59	0	0	0	0	0	0	0	0	1	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
60	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
61	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
62	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
63	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
64	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
65	0	0	0	6	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
66	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
67	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
68	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
69	0	0	1	0	0	0	0	0	0	Root;Bacteria
70	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
71	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
72	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
73	0	0	0	0	0	5	0	0	0	Root;Bacteria;Bacteroidetes
74	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
75	1	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
76	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
77	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
78	1	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
79	2	3	8	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
80	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
81	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
82	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
83	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
84	1	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Ruminococcus
85	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
86	0	0	0	0	0	0	0	1	0	Root;Bacteria
87	0	0	1	0	0	2	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
88	0	0	0	0	0	0	0	1	0	Root;Bacteria
89	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
90	0	0	0	9	0	0	3	0	0	Root;Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales;Erysipelotrichaceae;Turicibacter
91	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Butyrivibrio
92	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
93	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes
94	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
95	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes
96	0	0	0	1	0	1	0	1	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
97	0	0	0	0	0	1	0	0	0	Root;Bacteria
98	0	0	0	0	0	0	0	1	0	Root;Bacteria
99	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
100	0	0	0	1	0	0	0	0	0	Root;Bacteria
101	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
102	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
103	0	1	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
104	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
105	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
106	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
107	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
108	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Incertae Sedis XIII;Anaerovorax
109	0	0	0	1	0	0	1	5	2	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
110	0	0	0	0	0	2	0	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Olsenella
111	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
112	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
113	0	0	0	0	0	1	0	0	0	Root;Bacteria
114	0	0	0	0	0	1	0	0	0	Root;Bacteria
115	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
116	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
117	1	0	2	0	0	6	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
118	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
119	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
120	1	3	1	2	1	9	2	4	5	Root;Bacteria;Bacteroidetes
121	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
122	0	0	0	1	0	2	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
123	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
124	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
125	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
126	0	0	2	0	0	0	0	1	0	Root;Bacteria
127	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
128	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
129	0	0	0	1	0	0	0	0	0	Root;Bacteria
130	0	0	0	0	5	2	0	0	0	Root;Bacteria;Proteobacteria;Epsilonproteobacteria;Campylobacterales;Helicobacteraceae;Helicobacter
131	0	0	1	3	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
132	0	0	0	0	1	0	0	0	0	Root;Bacteria
133	0	0	1	0	0	0	0	0	0	Root;Bacteria
134	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
135	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
136	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
137	0	0	0	0	0	0	0	1	0	Root;Bacteria
138	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
139	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
140	0	0	0	0	0	0	1	3	0	Root;Bacteria
141	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
142	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
143	0	0	1	0	0	0	0	0	0	Root;Bacteria
144	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
145	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
146	1	0	0	0	2	0	2	0	3	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
147	0	1	0	1	1	0	0	0	3	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
148	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
149	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
150	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
151	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
152	0	0	0	1	0	0	1	2	19	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
153	0	2	1	2	0	0	1	1	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
154	2	18	0	1	0	0	21	4	4	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
155	0	0	0	0	0	5	9	5	3	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
156	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
157	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
158	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
159	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
160	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
161	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
162	0	0	0	0	0	3	5	2	6	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
163	0	0	0	0	0	0	0	0	1	Root;Bacteria
164	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
165	2	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
166	0	0	0	0	0	0	0	1	0	Root;Bacteria
167	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
168	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
169	0	2	0	7	0	0	0	2	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
170	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
171	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
172	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
173	0	0	0	0	0	1	0	0	0	Root;Bacteria
174	1	0	0	0	10	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Peptostreptococcaceae;Peptostreptococcaceae Incertae Sedis
175	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
176	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
177	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia
178	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
179	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
180	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
181	1	4	2	6	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
182	0	0	0	0	0	1	0	0	0	Root;Bacteria
183	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia
184	0	0	0	1	0	0	3	1	0	Root;Bacteria;Bacteroidetes
185	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
186	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
187	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
188	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
189	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
190	0	0	0	0	0	0	0	1	0	Root;Bacteria
191	2	1	10	2	24	0	0	1	1	Root;Bacteria
192	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus
193	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Butyrivibrio
194	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Acetanaerobacterium
195	0	0	0	0	0	1	0	0	0	Root;Bacteria
196	0	0	0	0	0	1	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
197	0	1	0	0	0	0	0	0	0	Root;Bacteria
198	0	2	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
199	0	0	0	0	0	1	1	0	0	Root;Bacteria
200	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
201	0	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
202	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
203	0	2	2	4	0	5	1	5	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
204	1	4	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
205	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
206	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
207	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
208	0	2	0	2	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
209	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
210	0	0	0	0	0	0	0	0	1	Root;Bacteria
211	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
212	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
213	0	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes
214	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
215	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
216	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
217	0	0	0	0	0	2	0	1	0	Root;Bacteria
218	0	0	0	0	9	1	0	0	0	Root;Bacteria;Bacteroidetes
219	0	0	0	0	1	0	0	0	0	Root;Bacteria
220	1	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
221	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes
222	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
223	0	0	0	0	0	0	0	2	2	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
224	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
225	0	2	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
226	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
227	0	1	2	0	9	1	1	1	3	Root;Bacteria;Bacteroidetes
228	16	0	0	0	12	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
229	0	0	0	0	0	1	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Incertae Sedis XIII
230	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
231	0	19	2	0	2	0	3	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
232	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
233	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
234	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus
235	0	1	1	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
236	0	0	0	0	0	2	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
237	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
238	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
239	0	0	0	0	0	1	0	0	0	Root;Bacteria
240	0	0	0	0	0	1	0	0	0	Root;Bacteria
241	0	0	0	0	0	0	2	0	0	Root;Bacteria;TM7;TM7_genera_incertae_sedis
242	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
243	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
244	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
245	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
246	0	0	0	0	0	0	0	1	0	Root;Bacteria
247	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
248	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus
249	1	0	0	0	0	0	0	0	0	Root;Bacteria
250	1	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
251	0	0	0	1	4	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
252	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
253	0	0	0	0	2	0	0	5	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
254	11	13	6	13	2	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
255	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
256	0	0	0	0	0	0	1	0	0	Root;Bacteria
257	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
258	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
259	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
260	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
261	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
262	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Bryantella
263	0	0	0	0	1	0	0	0	0	Root;Bacteria
264	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
265	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
266	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
267	1	0	0	5	17	20	0	0	0	Root;Bacteria
268	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
269	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
270	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
271	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
272	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
273	0	0	0	0	0	0	1	0	0	Root;Bacteria
274	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
275	0	0	0	0	0	0	1	0	0	Root;Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia
276	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
277	1	0	0	0	0	0	0	0	0	Root;Bacteria
278	0	0	0	0	0	1	0	0	0	Root;Bacteria
279	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
280	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
281	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
282	0	0	0	0	0	0	2	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
283	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
284	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
285	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
286	0	2	3	1	4	0	5	0	4	Root;Bacteria;Bacteroidetes
287	0	0	0	0	0	0	1	1	1	Root;Bacteria;Bacteroidetes
288	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
289	0	0	0	0	3	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
290	0	0	0	0	0	0	0	0	2	Root;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus
291	0	0	0	0	1	0	0	0	0	Root;Bacteria
292	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
293	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
294	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
295	29	1	10	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
296	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
297	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
298	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria
299	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
300	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia
301	0	0	0	0	0	0	2	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
302	0	0	0	0	0	1	0	0	0	Root;Bacteria
303	0	0	0	0	0	0	0	0	1	Root;Bacteria
304	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
305	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
306	0	0	0	0	0	0	0	0	1	Root;Bacteria
307	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
308	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Ruminococcaceae Incertae Sedis
309	0	0	0	1	0	0	0	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Denitrobacterium
310	0	0	1	0	0	0	0	0	0	Root;Bacteria
311	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
312	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
313	0	1	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
314	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
315	1	3	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
316	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
317	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
318	0	0	0	0	0	1	0	0	0	Root;Bacteria;Proteobacteria
319	0	2	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
320	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
321	0	0	0	0	0	0	0	0	1	Root;Bacteria
322	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
323	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
324	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
325	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
326	0	0	0	0	4	0	0	0	2	Root;Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales;Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
327	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
328	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
329	2	2	0	1	0	0	0	0	0	Root;Bacteria;Bacteroidetes
330	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes
331	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes
332	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
333	0	0	0	0	0	6	0	3	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
334	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
335	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
336	0	0	1	0	0	0	0	0	0	Root;Bacteria
337	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
338	0	0	0	0	0	0	0	1	0	Root;Bacteria
339	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
340	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
341	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
342	0	0	0	0	0	1	0	0	0	Root;Bacteria
343	0	0	0	0	0	0	0	0	1	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
344	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
345	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
346	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
347	0	0	0	1	0	0	0	0	0	Root;Bacteria
348	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
349	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
350	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
351	0	0	0	0	2	2	1	4	1	Root;Bacteria;Bacteroidetes
352	3	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
353	0	4	4	0	1	2	0	2	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
354	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
355	0	0	0	0	0	0	0	1	0	Root;Bacteria
356	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
357	0	0	0	4	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
358	0	0	1	0	0	0	0	0	0	Root;Bacteria
359	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
360	0	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
361	2	0	2	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
362	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
363	0	0	0	0	0	1	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae
364	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
365	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
366	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Roseburia
367	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
368	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
369	0	0	0	0	0	1	0	0	0	Root;Bacteria
370	2	1	0	5	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
371	1	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
372	0	1	0	0	0	0	0	0	0	Root;Bacteria
373	0	1	0	0	0	0	3	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridiaceae 1;Clostridium
374	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
375	0	0	0	0	0	0	4	0	0	Root;Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales;Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
376	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
377	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
378	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
379	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae
380	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus
381	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
382	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
383	4	9	0	2	0	0	0	2	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
384	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
385	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Carnobacteriaceae;Carnobacteriaceae 1
386	0	0	1	0	0	0	0	0	0	Root;Bacteria
387	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
388	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
389	0	1	0	0	0	0	0	0	0	Root;Bacteria
390	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
391	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes
392	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
393	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
394	0	0	1	0	0	0	0	0	0	Root;Bacteria
395	1	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
396	2	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
397	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
398	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
399	0	0	0	0	0	0	13	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
400	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
401	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
402	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
403	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Prevotellaceae
404	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae Incertae Sedis
405	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
406	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
407	1	0	0	0	0	4	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
408	1	5	3	2	0	0	0	0	1	Root;Bacteria;Bacteroidetes
409	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
410	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
411	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
412	0	0	0	0	2	0	0	0	0	Root;Bacteria;Bacteroidetes
413	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales
414	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae
415	0	0	0	0	0	7	0	2	2	Root;Bacteria;Bacteroidetes
416	0	1	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;Clostridia;Clostridiales"""

classic_otu_table1_no_tax = """#Full OTU Counts
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
0	0	0	0	0	0	0	0	1	0
1	0	0	0	0	0	1	0	0	0
2	0	0	0	0	0	0	0	0	1
3	2	1	0	0	0	0	0	0	0
4	1	0	0	0	0	0	0	0	0
5	0	0	0	0	0	0	0	0	1
6	0	0	0	0	0	0	0	1	0
7	0	0	2	0	0	0	0	0	2
8	1	1	0	2	4	0	0	0	0
9	0	0	2	0	0	0	0	0	0
10	0	1	0	0	0	0	0	0	0
11	0	0	0	0	0	0	1	0	0
12	0	0	0	0	0	0	1	0	0
13	1	0	0	1	0	1	0	0	0
14	0	0	1	1	0	0	0	0	0
15	0	0	0	0	1	0	0	0	0
16	1	0	2	0	0	0	0	0	0
17	0	0	0	1	0	0	4	10	37
18	0	1	0	0	0	0	0	0	0
19	0	0	0	0	0	0	0	0	1
20	0	0	0	0	1	0	0	0	0
21	0	0	0	0	0	0	2	3	2
22	0	0	0	0	2	0	1	0	0
23	14	1	14	1	0	0	0	0	0
24	1	0	0	0	0	0	0	0	0
25	0	0	0	1	0	0	0	0	0
26	0	0	0	0	0	0	0	1	1
27	0	0	0	0	0	0	0	0	1
28	0	1	0	0	0	0	0	0	0
29	6	0	4	0	2	0	0	0	0
30	0	0	0	0	0	1	0	0	0
31	1	0	0	0	0	0	0	0	0
32	0	0	0	0	1	0	0	0	0
33	0	0	0	1	0	0	0	0	0
34	0	0	0	0	0	0	8	10	2
35	1	0	1	0	0	0	0	0	0
36	1	0	1	0	0	0	0	1	1
37	0	0	0	0	0	1	0	0	0
38	0	0	1	0	0	0	0	0	0
39	0	0	0	0	0	0	0	1	0
40	0	0	1	0	0	0	0	0	0
41	0	0	1	0	0	0	0	1	0
42	0	0	0	0	0	1	0	0	0
43	0	0	0	0	0	1	0	0	0
44	0	0	1	0	0	0	0	0	0
45	1	0	0	0	0	0	0	0	0
46	0	0	0	0	0	0	0	0	1
47	0	0	0	1	0	0	0	0	0
48	0	0	0	0	1	0	0	0	0
49	0	0	0	1	0	0	0	0	0
50	0	1	0	0	0	0	0	0	0
51	0	1	0	0	0	0	0	0	0
52	0	2	0	0	0	0	0	0	0
53	0	0	0	0	0	0	2	0	1
54	0	0	0	0	0	0	5	0	0
55	0	0	0	0	0	0	1	0	0
56	0	0	0	0	0	1	0	0	0
57	0	0	0	0	0	0	0	1	0
58	1	0	1	0	0	0	0	0	0
59	0	0	0	0	0	0	0	0	1
60	0	0	0	0	0	0	0	1	0
61	0	0	1	0	0	0	0	1	0
62	0	0	1	0	0	0	0	0	0
63	1	0	1	0	0	0	0	0	0
64	0	0	0	0	0	0	0	0	1
65	0	0	0	6	0	0	0	1	0
66	0	0	1	0	0	0	0	0	0
67	0	0	1	0	0	0	0	0	0
68	1	0	0	0	0	0	0	0	0
69	0	0	1	0	0	0	0	0	0
70	0	0	0	0	0	1	0	0	0
71	0	0	1	0	0	0	0	0	0
72	0	0	0	0	0	1	0	0	0
73	0	0	0	0	0	5	0	0	0
74	0	0	0	1	0	0	0	0	0
75	1	0	1	0	0	0	0	0	0
76	0	0	0	1	0	0	0	0	0
77	0	0	0	1	0	0	0	0	0
78	1	0	1	1	0	0	0	0	0
79	2	3	8	0	1	0	0	0	0
80	0	0	0	0	0	0	0	0	1
81	1	0	0	0	0	0	0	0	0
82	0	0	0	0	0	2	0	0	0
83	0	0	0	1	0	0	0	1	0
84	1	0	0	0	0	0	0	2	0
85	0	0	0	0	0	0	0	0	1
86	0	0	0	0	0	0	0	1	0
87	0	0	1	0	0	2	0	1	0
88	0	0	0	0	0	0	0	1	0
89	0	0	1	0	0	0	0	0	0
90	0	0	0	9	0	0	3	0	0
91	0	0	0	1	0	0	0	0	0
92	0	0	0	0	0	0	1	0	0
93	0	0	0	0	0	0	2	1	0
94	0	0	0	0	0	0	0	1	0
95	0	0	0	2	0	0	0	0	0
96	0	0	0	1	0	1	0	1	1
97	0	0	0	0	0	1	0	0	0
98	0	0	0	0	0	0	0	1	0
99	0	0	0	1	0	0	0	0	0
100	0	0	0	1	0	0	0	0	0
101	0	0	0	3	0	0	0	0	0
102	0	1	0	0	0	0	0	0	0
103	0	1	0	0	0	0	1	0	0
104	0	0	0	0	0	1	0	0	0
105	0	1	0	0	0	0	0	0	0
106	0	0	0	0	0	1	0	0	0
107	0	0	0	0	0	1	0	0	0
108	0	0	0	0	0	0	1	0	0
109	0	0	0	1	0	0	1	5	2
110	0	0	0	0	0	2	0	0	0
111	0	0	0	0	0	0	1	0	0
112	0	0	0	0	0	0	1	0	0
113	0	0	0	0	0	1	0	0	0
114	0	0	0	0	0	1	0	0	0
115	0	0	0	0	0	1	0	0	0
116	0	1	0	0	0	0	0	0	0
117	1	0	2	0	0	6	0	0	0
118	0	0	0	1	0	0	0	0	0
119	0	0	0	0	0	0	0	1	0
120	1	3	1	2	1	9	2	4	5
121	0	0	0	0	0	0	0	1	0
122	0	0	0	1	0	2	0	0	0
123	0	0	0	0	0	0	1	0	0
124	0	0	0	0	0	0	1	0	0
125	0	0	0	0	0	0	1	0	0
126	0	0	2	0	0	0	0	1	0
127	0	0	0	0	0	1	0	0	0
128	0	0	0	0	0	0	1	0	0
129	0	0	0	1	0	0	0	0	0
130	0	0	0	0	5	2	0	0	0
131	0	0	1	3	0	0	0	0	0
132	0	0	0	0	1	0	0	0	0
133	0	0	1	0	0	0	0	0	0
134	0	0	0	0	0	0	0	0	1
135	0	0	1	0	0	0	0	0	0
136	1	0	0	0	0	0	0	0	0
137	0	0	0	0	0	0	0	1	0
138	0	0	1	0	0	0	0	0	0
139	1	0	0	0	0	0	0	0	0
140	0	0	0	0	0	0	1	3	0
141	0	0	0	0	1	0	0	0	0
142	0	0	0	0	1	0	0	0	0
143	0	0	1	0	0	0	0	0	0
144	0	0	0	0	0	1	0	0	0
145	0	0	2	0	0	0	0	0	0
146	1	0	0	0	2	0	2	0	3
147	0	1	0	1	1	0	0	0	3
148	0	0	0	0	0	1	0	0	0
149	0	0	0	0	0	0	0	1	0
150	0	0	0	0	1	0	0	0	0
151	0	0	0	1	0	0	0	1	0
152	0	0	0	1	0	0	1	2	19
153	0	2	1	2	0	0	1	1	1
154	2	18	0	1	0	0	21	4	4
155	0	0	0	0	0	5	9	5	3
156	0	0	1	0	0	0	0	1	0
157	0	0	1	0	0	0	0	0	0
158	1	0	1	0	0	0	0	0	0
159	0	0	0	0	0	0	0	1	1
160	0	0	0	0	0	0	1	0	0
161	0	0	1	0	0	0	0	0	0
162	0	0	0	0	0	3	5	2	6
163	0	0	0	0	0	0	0	0	1
164	0	0	0	0	0	1	0	0	0
165	2	1	1	0	0	0	0	0	0
166	0	0	0	0	0	0	0	1	0
167	1	0	0	0	0	0	0	0	0
168	0	0	0	1	0	0	0	0	0
169	0	2	0	7	0	0	0	2	0
170	0	0	0	1	0	0	0	0	0
171	0	0	0	1	0	0	0	0	0
172	1	0	0	0	0	0	0	0	0
173	0	0	0	0	0	1	0	0	0
174	1	0	0	0	10	0	0	0	0
175	0	0	0	0	1	0	0	0	0
176	0	0	0	0	0	1	0	0	0
177	0	0	0	1	0	0	0	0	0
178	0	0	0	2	0	0	0	0	0
179	0	0	0	1	0	0	0	0	0
180	0	0	0	0	1	0	0	0	0
181	1	4	2	6	0	0	0	0	0
182	0	0	0	0	0	1	0	0	0
183	0	0	0	0	0	0	1	0	0
184	0	0	0	1	0	0	3	1	0
185	0	0	0	0	0	0	0	0	1
186	0	0	1	0	0	0	0	0	0
187	0	1	0	0	0	0	0	0	1
188	0	0	0	0	0	0	0	1	0
189	0	0	0	1	0	0	0	0	0
190	0	0	0	0	0	0	0	1	0
191	2	1	10	2	24	0	0	1	1
192	0	0	0	0	0	1	0	0	0
193	0	0	0	0	0	1	0	0	0
194	0	0	2	0	0	0	0	0	0
195	0	0	0	0	0	1	0	0	0
196	0	0	0	0	0	1	0	1	0
197	0	1	0	0	0	0	0	0	0
198	0	2	0	0	0	1	0	0	0
199	0	0	0	0	0	1	1	0	0
200	0	0	0	2	0	0	0	0	0
201	0	0	0	1	0	1	0	0	0
202	0	0	0	0	0	0	1	0	0
203	0	2	2	4	0	5	1	5	0
204	1	4	0	1	0	0	0	0	0
205	0	0	0	0	0	0	0	1	0
206	0	1	0	0	0	0	0	0	0
207	0	0	0	0	0	0	0	1	0
208	0	2	0	2	0	0	0	1	0
209	0	0	1	0	0	0	0	0	0
210	0	0	0	0	0	0	0	0	1
211	1	0	0	1	0	0	0	0	0
212	0	0	0	0	0	0	0	0	1
213	0	0	0	0	0	0	0	2	0
214	0	0	0	0	0	0	0	1	0
215	0	0	0	0	0	0	0	1	0
216	0	0	0	0	0	0	0	1	0
217	0	0	0	0	0	2	0	1	0
218	0	0	0	0	9	1	0	0	0
219	0	0	0	0	1	0	0	0	0
220	1	0	0	0	1	0	0	0	0
221	0	0	0	0	0	0	0	1	0
222	0	1	0	0	0	0	0	0	0
223	0	0	0	0	0	0	0	2	2
224	0	0	0	1	0	0	0	0	0
225	0	2	1	0	0	0	0	0	0
226	0	0	0	0	0	1	0	0	0
227	0	1	2	0	9	1	1	1	3
228	16	0	0	0	12	0	0	0	0
229	0	0	0	0	0	1	1	0	0
230	0	0	0	1	0	0	0	0	0
231	0	19	2	0	2	0	3	0	0
232	0	0	0	0	0	0	1	0	0
233	0	0	0	0	1	0	0	0	0
234	0	0	0	0	1	0	0	0	0
235	0	1	1	0	1	0	0	0	0
236	0	0	0	0	0	2	0	0	0
237	0	0	0	0	1	0	0	0	0
238	0	0	0	0	0	0	0	1	0
239	0	0	0	0	0	1	0	0	0
240	0	0	0	0	0	1	0	0	0
241	0	0	0	0	0	0	2	0	0
242	0	0	0	0	0	0	1	0	0
243	0	0	0	0	0	0	1	0	0
244	0	0	0	0	0	0	0	0	1
245	0	0	0	1	0	0	0	1	0
246	0	0	0	0	0	0	0	1	0
247	0	0	1	0	0	0	0	0	0
248	1	0	0	1	0	0	0	0	0
249	1	0	0	0	0	0	0	0	0
250	1	0	0	0	0	0	0	1	0
251	0	0	0	1	4	0	0	0	0
252	0	0	0	1	0	0	0	0	0
253	0	0	0	0	2	0	0	5	0
254	11	13	6	13	2	0	0	0	0
255	0	0	0	0	0	1	0	0	0
256	0	0	0	0	0	0	1	0	0
257	0	0	0	0	0	0	5	0	0
258	0	0	1	0	0	0	0	0	0
259	0	0	0	0	0	0	0	1	0
260	0	0	0	0	0	0	0	1	0
261	0	0	0	0	0	0	0	1	0
262	0	1	0	0	0	0	0	0	1
263	0	0	0	0	1	0	0	0	0
264	0	0	0	0	0	1	0	0	0
265	0	0	0	0	0	2	0	0	0
266	0	0	0	2	0	0	0	0	0
267	1	0	0	5	17	20	0	0	0
268	0	0	0	0	0	0	1	0	0
269	0	0	0	1	0	0	0	0	0
270	0	0	1	0	0	0	0	0	0
271	0	0	0	0	0	0	0	0	1
272	0	0	0	1	0	0	0	0	0
273	0	0	0	0	0	0	1	0	0
274	0	0	0	0	0	0	1	0	0
275	0	0	0	0	0	0	1	0	0
276	0	0	0	0	0	0	0	1	0
277	1	0	0	0	0	0	0	0	0
278	0	0	0	0	0	1	0	0	0
279	0	0	0	0	0	1	0	0	0
280	0	1	0	0	0	0	0	0	0
281	1	0	0	0	0	0	0	0	0
282	0	0	0	0	0	0	2	0	0
283	0	0	0	0	0	0	2	1	0
284	0	0	0	1	0	0	0	0	0
285	0	0	0	0	0	0	1	0	0
286	0	2	3	1	4	0	5	0	4
287	0	0	0	0	0	0	1	1	1
288	0	0	0	0	0	1	0	0	0
289	0	0	0	0	3	0	0	0	0
290	0	0	0	0	0	0	0	0	2
291	0	0	0	0	1	0	0	0	0
292	0	0	0	0	1	0	0	0	0
293	0	0	0	0	0	1	0	0	0
294	0	1	0	0	0	0	0	0	0
295	29	1	10	0	0	0	0	0	0
296	0	0	0	0	1	0	0	0	0
297	0	0	0	1	0	0	0	0	0
298	0	0	0	0	0	0	1	0	0
299	0	0	0	0	0	0	1	0	1
300	0	0	0	0	0	1	0	0	0
301	0	0	0	0	0	0	2	0	0
302	0	0	0	0	0	1	0	0	0
303	0	0	0	0	0	0	0	0	1
304	0	0	0	0	0	0	0	1	0
305	1	0	0	0	0	0	0	0	0
306	0	0	0	0	0	0	0	0	1
307	0	0	1	0	0	0	0	0	0
308	0	1	0	0	0	0	0	0	0
309	0	0	0	1	0	0	0	0	0
310	0	0	1	0	0	0	0	0	0
311	0	0	0	0	0	1	0	0	0
312	0	0	1	0	0	0	0	0	0
313	0	1	0	0	0	0	0	0	1
314	0	0	1	0	0	0	0	0	0
315	1	3	1	0	0	0	0	0	0
316	0	1	0	0	0	0	0	0	0
317	0	0	0	0	0	0	1	0	0
318	0	0	0	0	0	1	0	0	0
319	0	2	1	0	0	0	0	0	0
320	0	0	0	1	0	0	0	0	0
321	0	0	0	0	0	0	0	0	1
322	0	0	0	1	0	0	0	0	0
323	0	0	1	0	0	0	0	0	0
324	0	0	1	0	0	0	0	0	0
325	0	1	0	0	0	0	0	0	0
326	0	0	0	0	4	0	0	0	2
327	0	0	0	0	0	0	0	1	0
328	0	0	0	1	0	0	0	0	0
329	2	2	0	1	0	0	0	0	0
330	0	0	1	0	0	0	0	0	0
331	0	0	0	0	1	0	0	0	0
332	0	1	0	0	0	0	0	0	0
333	0	0	0	0	0	6	0	3	0
334	1	0	0	0	0	0	0	0	0
335	0	0	0	0	0	0	0	1	0
336	0	0	1	0	0	0	0	0	0
337	0	0	0	1	0	0	0	0	0
338	0	0	0	0	0	0	0	1	0
339	0	0	1	0	0	0	0	0	0
340	0	0	2	0	0	0	0	0	0
341	0	0	1	0	0	0	0	1	0
342	0	0	0	0	0	1	0	0	0
343	0	0	0	0	0	0	0	0	1
344	0	0	1	0	0	0	0	0	0
345	1	0	0	0	0	0	0	0	0
346	0	1	0	0	0	0	0	0	0
347	0	0	0	1	0	0	0	0	0
348	0	0	0	1	0	0	0	0	0
349	0	0	0	0	0	0	1	0	1
350	1	0	0	0	0	0	0	0	0
351	0	0	0	0	2	2	1	4	1
352	3	0	0	0	0	0	0	0	0
353	0	4	4	0	1	2	0	2	1
354	0	0	0	0	0	1	0	0	0
355	0	0	0	0	0	0	0	1	0
356	0	0	0	0	0	1	0	0	0
357	0	0	0	4	0	0	0	0	0
358	0	0	1	0	0	0	0	0	0
359	0	0	1	0	0	0	0	0	0
360	0	0	1	0	0	0	0	1	1
361	2	0	2	1	0	0	0	0	0
362	1	0	0	1	0	0	0	0	0
363	0	0	0	0	0	1	0	1	0
364	1	0	0	0	0	0	0	0	0
365	0	0	0	0	0	2	0	0	0
366	0	0	0	1	0	0	0	0	0
367	0	0	0	0	1	0	0	0	0
368	0	0	0	0	0	1	0	0	0
369	0	0	0	0	0	1	0	0	0
370	2	1	0	5	0	1	0	0	0
371	1	1	0	0	0	0	0	0	0
372	0	1	0	0	0	0	0	0	0
373	0	1	0	0	0	0	3	0	0
374	0	0	0	0	0	0	1	0	0
375	0	0	0	0	0	0	4	0	0
376	0	0	0	0	0	0	0	0	1
377	0	0	0	0	0	0	0	1	0
378	0	0	0	0	0	0	0	0	1
379	0	0	0	0	0	1	0	0	0
380	0	0	0	0	0	0	0	0	1
381	0	0	2	0	0	0	0	0	0
382	0	0	0	0	0	0	0	1	0
383	4	9	0	2	0	0	0	2	0
384	0	0	1	0	0	0	0	0	0
385	0	0	0	0	0	0	0	0	1
386	0	0	1	0	0	0	0	0	0
387	0	0	1	0	0	0	0	0	0
388	0	0	0	1	0	0	0	0	0
389	0	1	0	0	0	0	0	0	0
390	0	0	0	0	0	0	0	0	1
391	0	0	0	0	0	0	0	0	1
392	0	1	0	0	0	0	0	0	0
393	0	0	0	0	0	1	0	0	0
394	0	0	1	0	0	0	0	0	0
395	1	1	1	0	0	0	0	0	0
396	2	0	0	0	0	0	0	0	0
397	0	0	0	0	0	0	0	1	0
398	0	0	0	0	0	0	0	1	0
399	0	0	0	0	0	0	13	0	0
400	0	0	0	0	0	0	1	0	0
401	0	1	0	0	0	0	0	0	0
402	0	1	0	0	0	0	0	0	0
403	0	0	0	0	0	0	0	1	0
404	0	0	0	0	0	0	0	1	0
405	0	0	0	0	0	0	0	1	0
406	0	0	0	0	0	1	0	0	0
407	1	0	0	0	0	4	0	0	0
408	1	5	3	2	0	0	0	0	1
409	0	0	0	0	0	0	0	1	1
410	0	0	0	0	1	0	0	0	0
411	0	0	0	1	0	0	0	0	0
412	0	0	0	0	2	0	0	0	0
413	0	0	0	0	0	0	0	1	0
414	1	0	1	0	0	0	0	0	0
415	0	0	0	0	0	7	0	2	2
416	0	1	0	0	1	0	0	0	0"""

biom_otu_table1_w_tax = """{"rows": [{"id": "0", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "1", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "2", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "3", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "4", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "5", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "6", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria"]}}, {"id": "7", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "8", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "9", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "10", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "11", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "12", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "13", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "14", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "15", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "16", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "17", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "18", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "19", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "20", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "21", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "22", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "23", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "24", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "25", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "26", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "27", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "28", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "29", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "30", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "31", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "32", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "33", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "34", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "35", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "36", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "37", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "38", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "39", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "40", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "41", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "42", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "43", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "44", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "45", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Coprobacillus"]}}, {"id": "46", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "47", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "48", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "49", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "50", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "51", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "52", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "53", "metadata": {"taxonomy": ["Root", "Bacteria", "Proteobacteria", "Deltaproteobacteria"]}}, {"id": "54", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "55", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "56", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "57", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "58", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "59", "metadata": {"taxonomy": ["Root", "Bacteria", "Deferribacteres", "Deferribacteres", "Deferribacterales", "Deferribacteraceae", "Mucispirillum"]}}, {"id": "60", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "61", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "62", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "63", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "64", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "65", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "66", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "67", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "68", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "69", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "70", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "71", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "72", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "73", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "74", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "75", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "76", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "77", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "78", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "79", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "80", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "81", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "82", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "83", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "84", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae", "Ruminococcus"]}}, {"id": "85", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "86", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "87", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "88", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "89", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "90", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Turicibacter"]}}, {"id": "91", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Butyrivibrio"]}}, {"id": "92", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "93", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "94", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "95", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "96", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "97", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "98", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "99", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "100", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "101", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "102", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "103", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "104", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "105", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "106", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "107", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "108", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Incertae Sedis XIII", "Anaerovorax"]}}, {"id": "109", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "110", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae", "Olsenella"]}}, {"id": "111", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "112", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "113", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "114", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "115", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "116", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "117", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "118", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "119", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "120", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "121", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "122", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "123", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]}}, {"id": "124", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]}}, {"id": "125", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "126", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "127", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "128", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "129", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "130", "metadata": {"taxonomy": ["Root", "Bacteria", "Proteobacteria", "Epsilonproteobacteria", "Campylobacterales", "Helicobacteraceae", "Helicobacter"]}}, {"id": "131", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "132", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "133", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "134", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "135", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "136", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "137", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "138", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "139", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "140", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "141", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "142", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "143", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "144", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "145", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "146", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "147", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "148", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "149", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "150", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "151", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "152", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "153", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "154", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "155", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "156", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "157", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "158", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "159", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "160", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "161", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "162", "metadata": {"taxonomy": ["Root", "Bacteria", "Deferribacteres", "Deferribacteres", "Deferribacterales", "Deferribacteraceae", "Mucispirillum"]}}, {"id": "163", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "164", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "165", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "166", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "167", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "168", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "169", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "170", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "171", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "172", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "173", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "174", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Peptostreptococcaceae", "Peptostreptococcaceae Incertae Sedis"]}}, {"id": "175", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "176", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "177", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia"]}}, {"id": "178", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "179", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "180", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "181", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "182", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "183", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia"]}}, {"id": "184", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "185", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "186", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "187", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "188", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "189", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "190", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "191", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "192", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Streptococcaceae", "Streptococcus"]}}, {"id": "193", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Butyrivibrio"]}}, {"id": "194", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae", "Acetanaerobacterium"]}}, {"id": "195", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "196", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "197", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "198", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales"]}}, {"id": "199", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "200", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "201", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "202", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "203", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "204", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "205", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "206", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "207", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "208", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "209", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "210", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "211", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "212", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "213", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "214", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "215", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "216", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "217", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "218", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "219", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "220", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "221", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "222", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "223", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "224", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "225", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "226", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "227", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "228", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "229", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Incertae Sedis XIII"]}}, {"id": "230", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "231", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "232", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "233", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "234", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "235", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "236", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales"]}}, {"id": "237", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "238", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "239", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "240", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "241", "metadata": {"taxonomy": ["Root", "Bacteria", "TM7", "TM7_genera_incertae_sedis"]}}, {"id": "242", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "243", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "244", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "245", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "246", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "247", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "248", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "249", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "250", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "251", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "252", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "253", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "254", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "255", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "256", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "257", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "258", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "259", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "260", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "261", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "262", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Bryantella"]}}, {"id": "263", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "264", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "265", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "266", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "267", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "268", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "269", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "270", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "271", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "272", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "273", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "274", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "275", "metadata": {"taxonomy": ["Root", "Bacteria", "Verrucomicrobia", "Verrucomicrobiae", "Verrucomicrobiales", "Verrucomicrobiaceae", "Akkermansia"]}}, {"id": "276", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "277", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "278", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "279", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "280", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "281", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "282", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "283", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "284", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "285", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "286", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "287", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "288", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "289", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "290", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus"]}}, {"id": "291", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "292", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "293", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "294", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "295", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "296", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "297", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "298", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria"]}}, {"id": "299", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "300", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia"]}}, {"id": "301", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "302", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "303", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "304", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "305", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "306", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "307", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "308", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae", "Ruminococcaceae Incertae Sedis"]}}, {"id": "309", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae", "Denitrobacterium"]}}, {"id": "310", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "311", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "312", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "313", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "314", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "315", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "316", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "317", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "318", "metadata": {"taxonomy": ["Root", "Bacteria", "Proteobacteria"]}}, {"id": "319", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "320", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "321", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "322", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "323", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "324", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "325", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "326", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Erysipelotrichaceae Incertae Sedis"]}}, {"id": "327", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "328", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "329", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "330", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "331", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "332", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "333", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "334", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "335", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "336", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "337", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "338", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "339", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "340", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "341", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "342", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "343", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]}}, {"id": "344", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "345", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "346", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "347", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "348", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "349", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "350", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "351", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "352", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "353", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "354", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "355", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "356", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "357", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "358", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "359", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "360", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "361", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "362", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "363", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae"]}}, {"id": "364", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "365", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "366", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Roseburia"]}}, {"id": "367", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "368", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "369", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "370", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "371", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "372", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "373", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Clostridiaceae", "Clostridiaceae 1", "Clostridium"]}}, {"id": "374", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "375", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Erysipelotrichaceae Incertae Sedis"]}}, {"id": "376", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "377", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "378", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "379", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "380", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus"]}}, {"id": "381", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "382", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "383", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "384", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "385", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Carnobacteriaceae", "Carnobacteriaceae 1"]}}, {"id": "386", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "387", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "388", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "389", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "390", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "391", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "392", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "393", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "394", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "395", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "396", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "397", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "398", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "399", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "400", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "401", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "402", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "403", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Prevotellaceae"]}}, {"id": "404", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "405", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "406", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "407", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "408", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "409", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "410", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "411", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "412", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "413", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "414", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "415", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "416", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [2, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 2, 0, 0, 0, 0, 0, 2], [1, 1, 0, 2, 4, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [1, 0, 0, 1, 0, 1, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [1, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 4, 10, 37], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 3, 2], [0, 0, 0, 0, 2, 0, 1, 0, 0], [14, 1, 14, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0, 0, 0, 0], [6, 0, 4, 0, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 8, 10, 2], [1, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 1], [0, 0, 0, 0, 0, 0, 5, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 6, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 5, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 1, 1, 0, 0, 0, 0, 0], [2, 3, 8, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 2, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 9, 0, 0, 3, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 2, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 1, 0, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 3, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 1, 5, 2], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 2, 0, 0, 6, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 3, 1, 2, 1, 9, 2, 4, 5], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 2, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 2, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 5, 2, 0, 0, 0], [0, 0, 1, 3, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 3, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 2, 0, 2, 0, 3], [0, 1, 0, 1, 1, 0, 0, 0, 3], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 1, 2, 19], [0, 2, 1, 2, 0, 0, 1, 1, 1], [2, 18, 0, 1, 0, 0, 21, 4, 4], [0, 0, 0, 0, 0, 5, 9, 5, 3], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 3, 5, 2, 6], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [2, 1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 2, 0, 7, 0, 0, 0, 2, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 10, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [1, 4, 2, 6, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 3, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [2, 1, 10, 2, 24, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 2, 2, 4, 0, 5, 1, 5, 0], [1, 4, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 2, 0, 2, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 2, 0, 1, 0], [0, 0, 0, 0, 9, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 2, 2], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 2, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 2, 0, 9, 1, 1, 1, 3], [16, 0, 0, 0, 12, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 19, 2, 0, 2, 0, 3, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 1, 1, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 4, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 2, 0, 0, 5, 0], [11, 13, 6, 13, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 5, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [1, 0, 0, 5, 17, 20, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 0, 2, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 2, 3, 1, 4, 0, 5, 0, 4], [0, 0, 0, 0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 3, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 2], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [29, 1, 10, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 3, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 2, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 4, 0, 0, 0, 2], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [2, 2, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 6, 0, 3, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 2, 2, 1, 4, 1], [3, 0, 0, 0, 0, 0, 0, 0, 0], [0, 4, 4, 0, 1, 2, 0, 2, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 4, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 1], [2, 0, 2, 1, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [2, 1, 0, 5, 0, 1, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 3, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 4, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [4, 9, 0, 2, 0, 0, 0, 2, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 1, 1, 0, 0, 0, 0, 0, 0], [2, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 13, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 0, 4, 0, 0, 0], [1, 5, 3, 2, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 7, 0, 2, 2], [0, 1, 0, 0, 1, 0, 0, 0, 0]], "columns": [{"id": "PC.354", "metadata": null}, {"id": "PC.355", "metadata": null}, {"id": "PC.356", "metadata": null}, {"id": "PC.481", "metadata": null}, {"id": "PC.593", "metadata": null}, {"id": "PC.607", "metadata": null}, {"id": "PC.634", "metadata": null}, {"id": "PC.635", "metadata": null}, {"id": "PC.636", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2604", "matrix_type": "dense", "shape": [417, 9], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-22T11:28:38.642265", "type": "OTU table", "id": null, "matrix_element_type": "int"}"""

biom_otu_table1_no_tax = """{"rows": [{"id": "0", "metadata": null}, {"id": "1", "metadata": null}, {"id": "2", "metadata": null}, {"id": "3", "metadata": null}, {"id": "4", "metadata": null}, {"id": "5", "metadata": null}, {"id": "6", "metadata": null}, {"id": "7", "metadata": null}, {"id": "8", "metadata": null}, {"id": "9", "metadata": null}, {"id": "10", "metadata": null}, {"id": "11", "metadata": null}, {"id": "12", "metadata": null}, {"id": "13", "metadata": null}, {"id": "14", "metadata": null}, {"id": "15", "metadata": null}, {"id": "16", "metadata": null}, {"id": "17", "metadata": null}, {"id": "18", "metadata": null}, {"id": "19", "metadata": null}, {"id": "20", "metadata": null}, {"id": "21", "metadata": null}, {"id": "22", "metadata": null}, {"id": "23", "metadata": null}, {"id": "24", "metadata": null}, {"id": "25", "metadata": null}, {"id": "26", "metadata": null}, {"id": "27", "metadata": null}, {"id": "28", "metadata": null}, {"id": "29", "metadata": null}, {"id": "30", "metadata": null}, {"id": "31", "metadata": null}, {"id": "32", "metadata": null}, {"id": "33", "metadata": null}, {"id": "34", "metadata": null}, {"id": "35", "metadata": null}, {"id": "36", "metadata": null}, {"id": "37", "metadata": null}, {"id": "38", "metadata": null}, {"id": "39", "metadata": null}, {"id": "40", "metadata": null}, {"id": "41", "metadata": null}, {"id": "42", "metadata": null}, {"id": "43", "metadata": null}, {"id": "44", "metadata": null}, {"id": "45", "metadata": null}, {"id": "46", "metadata": null}, {"id": "47", "metadata": null}, {"id": "48", "metadata": null}, {"id": "49", "metadata": null}, {"id": "50", "metadata": null}, {"id": "51", "metadata": null}, {"id": "52", "metadata": null}, {"id": "53", "metadata": null}, {"id": "54", "metadata": null}, {"id": "55", "metadata": null}, {"id": "56", "metadata": null}, {"id": "57", "metadata": null}, {"id": "58", "metadata": null}, {"id": "59", "metadata": null}, {"id": "60", "metadata": null}, {"id": "61", "metadata": null}, {"id": "62", "metadata": null}, {"id": "63", "metadata": null}, {"id": "64", "metadata": null}, {"id": "65", "metadata": null}, {"id": "66", "metadata": null}, {"id": "67", "metadata": null}, {"id": "68", "metadata": null}, {"id": "69", "metadata": null}, {"id": "70", "metadata": null}, {"id": "71", "metadata": null}, {"id": "72", "metadata": null}, {"id": "73", "metadata": null}, {"id": "74", "metadata": null}, {"id": "75", "metadata": null}, {"id": "76", "metadata": null}, {"id": "77", "metadata": null}, {"id": "78", "metadata": null}, {"id": "79", "metadata": null}, {"id": "80", "metadata": null}, {"id": "81", "metadata": null}, {"id": "82", "metadata": null}, {"id": "83", "metadata": null}, {"id": "84", "metadata": null}, {"id": "85", "metadata": null}, {"id": "86", "metadata": null}, {"id": "87", "metadata": null}, {"id": "88", "metadata": null}, {"id": "89", "metadata": null}, {"id": "90", "metadata": null}, {"id": "91", "metadata": null}, {"id": "92", "metadata": null}, {"id": "93", "metadata": null}, {"id": "94", "metadata": null}, {"id": "95", "metadata": null}, {"id": "96", "metadata": null}, {"id": "97", "metadata": null}, {"id": "98", "metadata": null}, {"id": "99", "metadata": null}, {"id": "100", "metadata": null}, {"id": "101", "metadata": null}, {"id": "102", "metadata": null}, {"id": "103", "metadata": null}, {"id": "104", "metadata": null}, {"id": "105", "metadata": null}, {"id": "106", "metadata": null}, {"id": "107", "metadata": null}, {"id": "108", "metadata": null}, {"id": "109", "metadata": null}, {"id": "110", "metadata": null}, {"id": "111", "metadata": null}, {"id": "112", "metadata": null}, {"id": "113", "metadata": null}, {"id": "114", "metadata": null}, {"id": "115", "metadata": null}, {"id": "116", "metadata": null}, {"id": "117", "metadata": null}, {"id": "118", "metadata": null}, {"id": "119", "metadata": null}, {"id": "120", "metadata": null}, {"id": "121", "metadata": null}, {"id": "122", "metadata": null}, {"id": "123", "metadata": null}, {"id": "124", "metadata": null}, {"id": "125", "metadata": null}, {"id": "126", "metadata": null}, {"id": "127", "metadata": null}, {"id": "128", "metadata": null}, {"id": "129", "metadata": null}, {"id": "130", "metadata": null}, {"id": "131", "metadata": null}, {"id": "132", "metadata": null}, {"id": "133", "metadata": null}, {"id": "134", "metadata": null}, {"id": "135", "metadata": null}, {"id": "136", "metadata": null}, {"id": "137", "metadata": null}, {"id": "138", "metadata": null}, {"id": "139", "metadata": null}, {"id": "140", "metadata": null}, {"id": "141", "metadata": null}, {"id": "142", "metadata": null}, {"id": "143", "metadata": null}, {"id": "144", "metadata": null}, {"id": "145", "metadata": null}, {"id": "146", "metadata": null}, {"id": "147", "metadata": null}, {"id": "148", "metadata": null}, {"id": "149", "metadata": null}, {"id": "150", "metadata": null}, {"id": "151", "metadata": null}, {"id": "152", "metadata": null}, {"id": "153", "metadata": null}, {"id": "154", "metadata": null}, {"id": "155", "metadata": null}, {"id": "156", "metadata": null}, {"id": "157", "metadata": null}, {"id": "158", "metadata": null}, {"id": "159", "metadata": null}, {"id": "160", "metadata": null}, {"id": "161", "metadata": null}, {"id": "162", "metadata": null}, {"id": "163", "metadata": null}, {"id": "164", "metadata": null}, {"id": "165", "metadata": null}, {"id": "166", "metadata": null}, {"id": "167", "metadata": null}, {"id": "168", "metadata": null}, {"id": "169", "metadata": null}, {"id": "170", "metadata": null}, {"id": "171", "metadata": null}, {"id": "172", "metadata": null}, {"id": "173", "metadata": null}, {"id": "174", "metadata": null}, {"id": "175", "metadata": null}, {"id": "176", "metadata": null}, {"id": "177", "metadata": null}, {"id": "178", "metadata": null}, {"id": "179", "metadata": null}, {"id": "180", "metadata": null}, {"id": "181", "metadata": null}, {"id": "182", "metadata": null}, {"id": "183", "metadata": null}, {"id": "184", "metadata": null}, {"id": "185", "metadata": null}, {"id": "186", "metadata": null}, {"id": "187", "metadata": null}, {"id": "188", "metadata": null}, {"id": "189", "metadata": null}, {"id": "190", "metadata": null}, {"id": "191", "metadata": null}, {"id": "192", "metadata": null}, {"id": "193", "metadata": null}, {"id": "194", "metadata": null}, {"id": "195", "metadata": null}, {"id": "196", "metadata": null}, {"id": "197", "metadata": null}, {"id": "198", "metadata": null}, {"id": "199", "metadata": null}, {"id": "200", "metadata": null}, {"id": "201", "metadata": null}, {"id": "202", "metadata": null}, {"id": "203", "metadata": null}, {"id": "204", "metadata": null}, {"id": "205", "metadata": null}, {"id": "206", "metadata": null}, {"id": "207", "metadata": null}, {"id": "208", "metadata": null}, {"id": "209", "metadata": null}, {"id": "210", "metadata": null}, {"id": "211", "metadata": null}, {"id": "212", "metadata": null}, {"id": "213", "metadata": null}, {"id": "214", "metadata": null}, {"id": "215", "metadata": null}, {"id": "216", "metadata": null}, {"id": "217", "metadata": null}, {"id": "218", "metadata": null}, {"id": "219", "metadata": null}, {"id": "220", "metadata": null}, {"id": "221", "metadata": null}, {"id": "222", "metadata": null}, {"id": "223", "metadata": null}, {"id": "224", "metadata": null}, {"id": "225", "metadata": null}, {"id": "226", "metadata": null}, {"id": "227", "metadata": null}, {"id": "228", "metadata": null}, {"id": "229", "metadata": null}, {"id": "230", "metadata": null}, {"id": "231", "metadata": null}, {"id": "232", "metadata": null}, {"id": "233", "metadata": null}, {"id": "234", "metadata": null}, {"id": "235", "metadata": null}, {"id": "236", "metadata": null}, {"id": "237", "metadata": null}, {"id": "238", "metadata": null}, {"id": "239", "metadata": null}, {"id": "240", "metadata": null}, {"id": "241", "metadata": null}, {"id": "242", "metadata": null}, {"id": "243", "metadata": null}, {"id": "244", "metadata": null}, {"id": "245", "metadata": null}, {"id": "246", "metadata": null}, {"id": "247", "metadata": null}, {"id": "248", "metadata": null}, {"id": "249", "metadata": null}, {"id": "250", "metadata": null}, {"id": "251", "metadata": null}, {"id": "252", "metadata": null}, {"id": "253", "metadata": null}, {"id": "254", "metadata": null}, {"id": "255", "metadata": null}, {"id": "256", "metadata": null}, {"id": "257", "metadata": null}, {"id": "258", "metadata": null}, {"id": "259", "metadata": null}, {"id": "260", "metadata": null}, {"id": "261", "metadata": null}, {"id": "262", "metadata": null}, {"id": "263", "metadata": null}, {"id": "264", "metadata": null}, {"id": "265", "metadata": null}, {"id": "266", "metadata": null}, {"id": "267", "metadata": null}, {"id": "268", "metadata": null}, {"id": "269", "metadata": null}, {"id": "270", "metadata": null}, {"id": "271", "metadata": null}, {"id": "272", "metadata": null}, {"id": "273", "metadata": null}, {"id": "274", "metadata": null}, {"id": "275", "metadata": null}, {"id": "276", "metadata": null}, {"id": "277", "metadata": null}, {"id": "278", "metadata": null}, {"id": "279", "metadata": null}, {"id": "280", "metadata": null}, {"id": "281", "metadata": null}, {"id": "282", "metadata": null}, {"id": "283", "metadata": null}, {"id": "284", "metadata": null}, {"id": "285", "metadata": null}, {"id": "286", "metadata": null}, {"id": "287", "metadata": null}, {"id": "288", "metadata": null}, {"id": "289", "metadata": null}, {"id": "290", "metadata": null}, {"id": "291", "metadata": null}, {"id": "292", "metadata": null}, {"id": "293", "metadata": null}, {"id": "294", "metadata": null}, {"id": "295", "metadata": null}, {"id": "296", "metadata": null}, {"id": "297", "metadata": null}, {"id": "298", "metadata": null}, {"id": "299", "metadata": null}, {"id": "300", "metadata": null}, {"id": "301", "metadata": null}, {"id": "302", "metadata": null}, {"id": "303", "metadata": null}, {"id": "304", "metadata": null}, {"id": "305", "metadata": null}, {"id": "306", "metadata": null}, {"id": "307", "metadata": null}, {"id": "308", "metadata": null}, {"id": "309", "metadata": null}, {"id": "310", "metadata": null}, {"id": "311", "metadata": null}, {"id": "312", "metadata": null}, {"id": "313", "metadata": null}, {"id": "314", "metadata": null}, {"id": "315", "metadata": null}, {"id": "316", "metadata": null}, {"id": "317", "metadata": null}, {"id": "318", "metadata": null}, {"id": "319", "metadata": null}, {"id": "320", "metadata": null}, {"id": "321", "metadata": null}, {"id": "322", "metadata": null}, {"id": "323", "metadata": null}, {"id": "324", "metadata": null}, {"id": "325", "metadata": null}, {"id": "326", "metadata": null}, {"id": "327", "metadata": null}, {"id": "328", "metadata": null}, {"id": "329", "metadata": null}, {"id": "330", "metadata": null}, {"id": "331", "metadata": null}, {"id": "332", "metadata": null}, {"id": "333", "metadata": null}, {"id": "334", "metadata": null}, {"id": "335", "metadata": null}, {"id": "336", "metadata": null}, {"id": "337", "metadata": null}, {"id": "338", "metadata": null}, {"id": "339", "metadata": null}, {"id": "340", "metadata": null}, {"id": "341", "metadata": null}, {"id": "342", "metadata": null}, {"id": "343", "metadata": null}, {"id": "344", "metadata": null}, {"id": "345", "metadata": null}, {"id": "346", "metadata": null}, {"id": "347", "metadata": null}, {"id": "348", "metadata": null}, {"id": "349", "metadata": null}, {"id": "350", "metadata": null}, {"id": "351", "metadata": null}, {"id": "352", "metadata": null}, {"id": "353", "metadata": null}, {"id": "354", "metadata": null}, {"id": "355", "metadata": null}, {"id": "356", "metadata": null}, {"id": "357", "metadata": null}, {"id": "358", "metadata": null}, {"id": "359", "metadata": null}, {"id": "360", "metadata": null}, {"id": "361", "metadata": null}, {"id": "362", "metadata": null}, {"id": "363", "metadata": null}, {"id": "364", "metadata": null}, {"id": "365", "metadata": null}, {"id": "366", "metadata": null}, {"id": "367", "metadata": null}, {"id": "368", "metadata": null}, {"id": "369", "metadata": null}, {"id": "370", "metadata": null}, {"id": "371", "metadata": null}, {"id": "372", "metadata": null}, {"id": "373", "metadata": null}, {"id": "374", "metadata": null}, {"id": "375", "metadata": null}, {"id": "376", "metadata": null}, {"id": "377", "metadata": null}, {"id": "378", "metadata": null}, {"id": "379", "metadata": null}, {"id": "380", "metadata": null}, {"id": "381", "metadata": null}, {"id": "382", "metadata": null}, {"id": "383", "metadata": null}, {"id": "384", "metadata": null}, {"id": "385", "metadata": null}, {"id": "386", "metadata": null}, {"id": "387", "metadata": null}, {"id": "388", "metadata": null}, {"id": "389", "metadata": null}, {"id": "390", "metadata": null}, {"id": "391", "metadata": null}, {"id": "392", "metadata": null}, {"id": "393", "metadata": null}, {"id": "394", "metadata": null}, {"id": "395", "metadata": null}, {"id": "396", "metadata": null}, {"id": "397", "metadata": null}, {"id": "398", "metadata": null}, {"id": "399", "metadata": null}, {"id": "400", "metadata": null}, {"id": "401", "metadata": null}, {"id": "402", "metadata": null}, {"id": "403", "metadata": null}, {"id": "404", "metadata": null}, {"id": "405", "metadata": null}, {"id": "406", "metadata": null}, {"id": "407", "metadata": null}, {"id": "408", "metadata": null}, {"id": "409", "metadata": null}, {"id": "410", "metadata": null}, {"id": "411", "metadata": null}, {"id": "412", "metadata": null}, {"id": "413", "metadata": null}, {"id": "414", "metadata": null}, {"id": "415", "metadata": null}, {"id": "416", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [2, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 2, 0, 0, 0, 0, 0, 2], [1, 1, 0, 2, 4, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [1, 0, 0, 1, 0, 1, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [1, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 4, 10, 37], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 3, 2], [0, 0, 0, 0, 2, 0, 1, 0, 0], [14, 1, 14, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0, 0, 0, 0], [6, 0, 4, 0, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 8, 10, 2], [1, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 1], [0, 0, 0, 0, 0, 0, 5, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 6, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 5, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 1, 1, 0, 0, 0, 0, 0], [2, 3, 8, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 2, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 9, 0, 0, 3, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 2, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 1, 0, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 3, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 1, 5, 2], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 2, 0, 0, 6, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 3, 1, 2, 1, 9, 2, 4, 5], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 2, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 2, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 5, 2, 0, 0, 0], [0, 0, 1, 3, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 3, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 2, 0, 2, 0, 3], [0, 1, 0, 1, 1, 0, 0, 0, 3], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 1, 2, 19], [0, 2, 1, 2, 0, 0, 1, 1, 1], [2, 18, 0, 1, 0, 0, 21, 4, 4], [0, 0, 0, 0, 0, 5, 9, 5, 3], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 3, 5, 2, 6], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [2, 1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 2, 0, 7, 0, 0, 0, 2, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 10, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [1, 4, 2, 6, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 3, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [2, 1, 10, 2, 24, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 2, 2, 4, 0, 5, 1, 5, 0], [1, 4, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 2, 0, 2, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 2, 0, 1, 0], [0, 0, 0, 0, 9, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 2, 2], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 2, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 2, 0, 9, 1, 1, 1, 3], [16, 0, 0, 0, 12, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 19, 2, 0, 2, 0, 3, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 1, 1, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 4, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 2, 0, 0, 5, 0], [11, 13, 6, 13, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 5, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [1, 0, 0, 5, 17, 20, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 0, 2, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 2, 3, 1, 4, 0, 5, 0, 4], [0, 0, 0, 0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 3, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 2], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [29, 1, 10, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 3, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 2, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 4, 0, 0, 0, 2], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [2, 2, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 6, 0, 3, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 2, 2, 1, 4, 1], [3, 0, 0, 0, 0, 0, 0, 0, 0], [0, 4, 4, 0, 1, 2, 0, 2, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 4, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 1], [2, 0, 2, 1, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [2, 1, 0, 5, 0, 1, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 3, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 4, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [4, 9, 0, 2, 0, 0, 0, 2, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 1, 1, 0, 0, 0, 0, 0, 0], [2, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 13, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 0, 4, 0, 0, 0], [1, 5, 3, 2, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 7, 0, 2, 2], [0, 1, 0, 0, 1, 0, 0, 0, 0]], "columns": [{"id": "PC.354", "metadata": null}, {"id": "PC.355", "metadata": null}, {"id": "PC.356", "metadata": null}, {"id": "PC.481", "metadata": null}, {"id": "PC.593", "metadata": null}, {"id": "PC.607", "metadata": null}, {"id": "PC.634", "metadata": null}, {"id": "PC.635", "metadata": null}, {"id": "PC.636", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2604", "matrix_type": "dense", "shape": [417, 9], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-22T11:28:47.373680", "type": "OTU table", "id": null, "matrix_element_type": "int"}"""



if __name__ =='__main__':
    main()
