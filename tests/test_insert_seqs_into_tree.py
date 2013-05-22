#!/usr/bin/env python
# File created on 11 Oct 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from qiime.insert_seqs_into_tree import convert_tree_tips, \
                                    write_updated_tree_file, \
                                    strip_and_rename_unwanted_labels_from_tree
from os.path import splitext
from os import getcwd, remove, rmdir, mkdir
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from StringIO import StringIO
from qiime.util import get_tmp_filename

class Tests(TestCase):

    def setUp(self):
        '''setup the files for testing pplacer'''
        
        # create a list of files to cleanup
        self._paths_to_clean_up = []
        self._dirs_to_clean_up = []
        
        # get a tmp filename to use
        self.basename=splitext(get_tmp_filename())[0]
        
        self.align_map={'seq0000005': 'Species005', 'seq0000004': 'Species004',\
                        'seq0000007': 'Species007', 'seq0000006': 'Species006',\
                        'seq0000001': 'Species001', 'seq0000003': 'Species003',\
                        'seq0000002': 'Species002'}
        
        # create and write out RAxML stats file
        self.tmp_tree_fname=self.basename+'.tre'
        tree_out=open(self.tmp_tree_fname,'w')
        tree_out.write(STARTING_TREE)
        tree_out.close()
        self._paths_to_clean_up.append(self.tmp_tree_fname)

    def tearDown(self): 
        """cleans up all files initially created"""
        # remove the tempdir and contents
        map(remove,self._paths_to_clean_up)
        map(rmdir,self._dirs_to_clean_up)
    
class insertSeqsTests(Tests):
    """Tests for the pplacer application controller"""
    
    def test_convert_tree_tips(self):
        """Convert tree tips to phylip labels"""
        
        # convert tree tips to PHYLIP labels
        tree=convert_tree_tips(self.align_map,self.tmp_tree_fname)
        
        self.assertEqual(tree.getNewick(with_distances=True), \
                         PHYLIP_TREE)

        
    def test_write_updated_tree_file(self):
        """Write tree out"""
        
        # create temp filename
        new_tree_fp=splitext(get_tmp_filename())[0]+'.tre'
        self._paths_to_clean_up.append(new_tree_fp)
        
        # parse and load tree
        tree=DndParser(StringIO(STARTING_TREE), constructor=PhyloNode)
        
        # write out temp tree
        write_updated_tree_file(new_tree_fp,tree)
        
        self.assertTrue(open(new_tree_fp).read()>0)
        
    def test_strip_and_rename_unwanted_labels_from_tree(self):
        """Remove unwanted text from Tip labels"""
        
        # parse and load tree
        result=DndParser(StringIO(RESULTING_QUERY_TREE), constructor=PhyloNode)
        
        # strip and rename tips
        result_tree=strip_and_rename_unwanted_labels_from_tree(self.align_map,\
                                                               result)
        self.assertEqual(result_tree.getNewick(with_distances=True), \
                         STRIPPED_TREE)

STARTING_TREE="""(Species002:0.00000043418318065054,((Species003:0.01932550067944402081,Species004:0.08910446960529855298):0.00000043418318065054,Species005:0.17394765077611337722):0.00000043418318065054,Species001:0.00000043418318065054):0.0;"""

PHYLIP_TREE="""(seq0000002:4.34183180651e-07,((seq0000003:0.0193255006794,seq0000004:0.0891044696053):4.34183180651e-07,seq0000005:0.173947650776):4.34183180651e-07,seq0000001:4.34183180651e-07):0.0;"""

RESULTING_QUERY_TREE="""(seq0000003:1.0,(QUERY___seq0000006___1,QUERY___seq0000007___1,seq0000004:1.0):1.0,(seq0000005:1.0,(seq0000001:1.0,seq0000002:1.0):1.0):1.0);"""

STRIPPED_TREE="""(Species003:1.0,(Species006,Species007,Species004:1.0):1.0,(Species005:1.0,(Species001:1.0,Species002:1.0):1.0):1.0);"""

if __name__ == "__main__":
    main()