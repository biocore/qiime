#!/usr/bin/env python
#file test_filter_otus_by_sample.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"

from numpy import array
from StringIO import StringIO
from os.path import exists
from cogent.util.unit_test import TestCase, main
from os import remove
from cogent import LoadSeqs
import shutil
from qiime.filter_otus_by_sample import (filter_otus,filter_aln_by_otus,\
                                 process_extract_samples)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        self.aln=[]
        self.aln.append(('SampleA','AAAAAAAAAAAAAAA'))
        self.aln.append(('SampleB','CCCCCCC'))
        self.aln.append(('SampleC','GGGGGGGGGGGGGG'))
        
        self.otus={'0':['SampleC'],'1':['SampleC','SampleA'],'2':['SampleB',\
                                                                  'SampleA']}

        self.prefs={}
        self.prefs={'0':'SampleC','1':'SampleB'}

    def test_process_extract_samples(self):
        """process_extract_samples: parses the cmd line and determines which
samples should be removed"""

        self.sample_to_extract='SampleA,SampleB'
        exp1={'0':'SampleA','1':'SampleB'}

        obs1=process_extract_samples(self.sample_to_extract)

        self.assertEqual(obs1,exp1)

    def test_filter_otus(self):
        """filter_otus: determines which sequences should be removed and
generates a new otus list"""

        exp1=[('1',['SampleA']),('2',['SampleA'])]
        obs1=filter_otus(self.otus,self.prefs)
            
        self.assertEqual(obs1,exp1)

        
    def test_filter_aln_by_otus(self):
        """filter_aln_by_otus: determines which sequences to keep and which
sequences to remove"""

        self.sample_to_extract='SampleA,SampleB'
        exp1=[]
        exp1.append(('SampleA','AAAAAAAAAAAAAAA'))
        exp2=[]
        exp2.append(('SampleB','CCCCCCC'))
        exp2.append(('SampleC','GGGGGGGGGGGGGG'))
        aln=LoadSeqs(data=self.aln,aligned=False)

        obs1,obs2=filter_aln_by_otus(aln,self.prefs)

        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,exp2)
        

#run tests if called from command line
if __name__ == "__main__":
    main()
