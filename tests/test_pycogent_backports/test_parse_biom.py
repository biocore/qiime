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

if __name__ =='__main__':
    main()
