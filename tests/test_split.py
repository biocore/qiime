#!/usr/bin/env python
# File created on 24 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from qiime.split import (split_mapping_file_on_field, 
                         split_otu_table_on_sample_metadata)
from qiime.format import format_biom_table
from biom.parse import parse_biom_table
from biom.table import DenseOTUTable

class SplitTests(TestCase):
    """ Tests of the split module """

    def setUp(self):
        """ """
        self.mapping_f1 = mapping_f1.split('\n')
        self.mapping_f2 = mapping_f2.split('\n')
        self.mapping_exp = list(mapping_exp)
        self.otu_table_f1 = otu_table_f1.split('\n')
    
    def test_split_mapping_file_on_field(self):
        """ split_mapping_file_on_field functions as expected with valid input 
        """
        actual = list(split_mapping_file_on_field(self.mapping_f1,'Treatment'))
        actual.sort()
        self.mapping_exp.sort()
        self.assertEqual(actual,self.mapping_exp)
    
    def test_split_otu_table_on_sample_metadata(self):
        """ split_otu_table_on_sample_metadata functions as expected with valid input """
        actual = list(split_otu_table_on_sample_metadata(self.otu_table_f1,
                                                         self.mapping_f1,
                                                         "Treatment"))
        for id_, e in actual:
            try:
                parse_biom_table(e)
            except:
                print e
        actual = [(id_,parse_biom_table(e)) for id_, e in actual]
        exp = [(id_,parse_biom_table(e)) for id_, e in otu_table_exp1]
        
        actual.sort()
        exp.sort()
        
        for a,e in zip(actual,exp):
            self.assertEqual(a,e,"OTU tables are not equal:\n%s\n%s" % \
             (format_biom_table(a[1]),format_biom_table(e[1])))

    def test_split_otu_table_on_sample_metadata_extra_mapping_entries(self):
        """ split_otu_table_on_sample_metadata functions as expected with extra mapping data """
        actual = list(split_otu_table_on_sample_metadata(self.otu_table_f1,
                                                         self.mapping_f2,
                                                         "Treatment"))
                                                         
        actual = [(id_,parse_biom_table(e)) for id_, e in actual]
        exp = [(id_,parse_biom_table(e)) for id_, e in otu_table_exp1]
        
        actual.sort()
        exp.sort()
        
        for a,e in zip(actual,exp):
            self.assertEqual(a,e,"OTU tables are not equal:\n%s\n%s" % \
             (format_biom_table(a[1]),format_biom_table(e[1])))
        

mapping_f1 = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Co_ntrol	20061218	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Co_ntrol	20061126	Control_mouse_I.D._356
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	Control_mouse_I.D._481
PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	Control	20071210	Control_mouse_I.D._593
PC.607	AACTGTGCGTAC	YATGCTGCCTCCCGTAGGAGT	Fast	20071112	Fasting_mouse_I.D._607
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._636
"""

mapping_f2 = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Co_ntrol	20061218	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Co_ntrol	20061126	Control_mouse_I.D._356
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	Control_mouse_I.D._481
PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	Control	20071210	Control_mouse_I.D._593
PC.607	AACTGTGCGTAC	YATGCTGCCTCCCGTAGGAGT	Fast	20071112	Fasting_mouse_I.D._607
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._636
Fake.sample	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Other	20080116	Fasting_mouse_I.D._636
"""

mapping_exp = [("Co_ntrol", """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Co_ntrol	20061218	Control_mouse_I.D._354
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Co_ntrol	20061126	Control_mouse_I.D._356"""),

("Control","""#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._355
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	Control_mouse_I.D._481
PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	Control	20071210	Control_mouse_I.D._593"""),

("Fast","""#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
PC.607	AACTGTGCGTAC	YATGCTGCCTCCCGTAGGAGT	Fast	20071112	Fasting_mouse_I.D._607
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._636""")]

otu_table_f1 = """{
    "id":null,
    "format": "Biological Observation Matrix v0.9",
    "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
    "type": "OTU table",
    "generated_by": "QIIME",
    "date": "2011-12-19T19:00:00",
    "rows":[
            {"id":"GG_OTU_1", "metadata":null},
            {"id":"GG_OTU_2", "metadata":null},
            {"id":"GG_OTU_3", "metadata":null},
            {"id":"GG_OTU_4", "metadata":null},
            {"id":"GG_OTU_5", "metadata":null}
        ],
    "columns": [
            {"id":"PC.354", "metadata":null},
            {"id":"PC.355", "metadata":null},
            {"id":"PC.481", "metadata":null},
            {"id":"PC.607", "metadata":null},
            {"id":"PC.635", "metadata":null},
            {"id":"PC.356", "metadata":null},
            {"id":"PC.636", "metadata":null}
        ],
    "matrix_type": "sparse",
    "matrix_element_type": "int",
    "shape": [5, 7],
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
            [4,2,1],
            [4,6,42]
           ]
}"""

otu_table_exp1 = [("Co_ntrol","""{
    "id":null,
    "format": "Biological Observation Matrix v0.9",
    "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
    "type": "OTU table",
    "generated_by": "QIIME",
    "date": "2011-12-19T19:00:00",
    "rows":[
            {"id":"GG_OTU_1", "metadata":null},
            {"id":"GG_OTU_2", "metadata":null},
            {"id":"GG_OTU_3", "metadata":null},
            {"id":"GG_OTU_4", "metadata":null},
            {"id":"GG_OTU_5", "metadata":null}
        ],
    "columns": [
            {"id":"PC.354", "metadata":null},
            {"id":"PC.356", "metadata":null}
        ],
    "matrix_type": "sparse",
    "matrix_element_type": "float",
    "shape": [5, 2],
    "data":[
            [1,0,5.0],
            [1,1,1.0],
            [3,0,2.0],
            [3,1,1.0]
           ]
}"""),
("Control","""{
    "id":null,
    "format": "Biological Observation Matrix v0.9",
    "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
    "type": "OTU table",
    "generated_by": "QIIME",
    "date": "2011-12-19T19:00:00",
    "rows":[
            {"id":"GG_OTU_1", "metadata":null},
            {"id":"GG_OTU_2", "metadata":null},
            {"id":"GG_OTU_3", "metadata":null},
            {"id":"GG_OTU_4", "metadata":null},
            {"id":"GG_OTU_5", "metadata":null}
        ],
    "columns": [
            {"id":"PC.355", "metadata":null},
            {"id":"PC.481", "metadata":null}
        ],
    "matrix_type": "sparse",
    "matrix_element_type": "float",
    "shape": [5, 2],
    "data":[[0,1,1.0],
            [1,0,1.0],
            [2,1,1.0],
            [3,0,1.0],
            [3,1,1.0],
            [4,0,1.0],
            [4,1,1.0]
           ]
}
"""),
("Fast","""{
    "id":null,
    "format": "Biological Observation Matrix v0.9",
    "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
    "type": "OTU table",
    "generated_by": "QIIME",
    "date": "2011-12-19T19:00:00",
    "rows":[
            {"id":"GG_OTU_1", "metadata":null},
            {"id":"GG_OTU_2", "metadata":null},
            {"id":"GG_OTU_3", "metadata":null},
            {"id":"GG_OTU_4", "metadata":null},
            {"id":"GG_OTU_5", "metadata":null}
        ],
    "columns": [
            {"id":"PC.607", "metadata":null},
            {"id":"PC.635", "metadata":null},
            {"id":"PC.636", "metadata":null}
        ],
    "matrix_type": "sparse",
    "matrix_element_type": "float",
    "shape": [5, 3],
    "data":[[1,0,2.0],
            [1,1,3.0],
            [2,0,4.0],
            [2,1,2.0],
            [4,2,42]
           ]
}""")]
if __name__ == "__main__":
    main()