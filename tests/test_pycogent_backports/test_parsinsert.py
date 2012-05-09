#!/usr/bin/env python
"""Tests for ParsInsert v1.03 application controller."""

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"

from shutil import rmtree
from cogent.util.unit_test import TestCase, main
from qiime.pycogent_backports.parsinsert import ParsInsert, \
                                                insert_sequences_into_tree
from cogent.core.alignment import Alignment
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.tree import DndParser
from cogent.core.moltype import DNA
from cogent.app.util import get_tmp_filename
from os.path import splitext
from os import getcwd, remove, rmdir, mkdir

class ParsInsertTests(TestCase):
    def setUp(self):
        
        # create a list of files to cleanup
        self._paths_to_clean_up = []
        self._dirs_to_clean_up = []
        
        # load query seqs
        self.seqs = Alignment(MinimalFastaParser(QUERY_SEQS.split()))
        
        # generate temp filename
        tmp_dir='/tmp'
        self.outfile = get_tmp_filename(tmp_dir)
        
        # create and write out reference sequence file
        self.outfasta=splitext(self.outfile)[0]+'.fasta'
        fastaout=open(self.outfasta,'w')
        fastaout.write(REF_SEQS)
        fastaout.close()
        self._paths_to_clean_up.append(self.outfasta)
        
        # create and write out starting tree file
        self.outtree=splitext(self.outfile)[0]+'.tree'
        treeout=open(self.outtree,'w')
        treeout.write(REF_TREE)
        treeout.close()
        self._paths_to_clean_up.append(self.outtree)
    
    def tearDown(self): 
        """cleans up all files initially created"""
        # remove the tempdir and contents
        map(remove,self._paths_to_clean_up)
        map(rmdir,self._dirs_to_clean_up)
    
    def test_base_command(self):
        """Base command-calls"""
        
        app = ParsInsert()
        self.assertEqual(app.BaseCommand, \
                         ''.join(['cd "',getcwd(),'/"; ','ParsInsert']))
        
    def test_change_working_dir(self):
        """Change working dir"""
        
        app = ParsInsert(WorkingDir='/tmp/ParsInsertTest')
        self.assertEqual(app.BaseCommand, \
                       ''.join(['cd "','/tmp/ParsInsertTest',\
                                '/"; ','ParsInsert']))
                                
        rmtree('/tmp/ParsInsertTest')

    def test_insert_sequences_into_tree(self):
        """Inserts sequences into Tree"""
        
        # define log fp
        log_fp='/tmp/parsinsert.log'
        self._paths_to_clean_up.append(log_fp)
        
        # define tax assignment values fp
        tax_assign_fp='/tmp/tax_assignments.log'
        self._paths_to_clean_up.append(tax_assign_fp)
        
        # set the reference alignment and starting tree
        param={
                '-t':self.outtree,
                '-s':self.outfasta,
                '-l':log_fp,
                '-o':tax_assign_fp
              }
        
        seqs, align_map = self.seqs.toPhylip()
        
        # insert sequences into tree
        tree = insert_sequences_into_tree(seqs, DNA, params=param)

        # rename tips back to query names
        for node in tree.tips():
            if node.Name in align_map:
                node.Name = align_map[node.Name]
                
        self.assertEqual(tree.getNewick(with_distances=True),exp_tree)


        
QUERY_SEQS= """\
>6
TGCATGTCAGTATAGCTTTGGTGAAACTGCGAATGGCTCATTAAATCAGT
>7
TGCATGTCAGTATAACTTTGGTGAAACTGCGAATGGCTCATTAAATCAGT
""" 

REF_SEQS= """\
>seq0000011
TGCATGTCAGTATAGCTTTAGTGAAACTGCGAATGGCTCATTAAATCAGT
>seq0000012
TGCATGTCAGTATAGCTTTAGTGAAACTGCGAATGGCTNNTTAAATCAGT
>seq0000013
TGCATGTCAGTATAGCATTAGTGAAACTGCGAATGGCTCATTAAATCAGT
>seq0000014
TCCATGTCAGTATAACTTTGGTGAAACTGCGAATGGCTCATTAAATCAGG
>seq0000015
NNNNNNNNNNTATATCTTATGTGAAACTTCGAATGCCTCATTAAATCAGT
"""

REF_TREE="""((seq0000014:0.08408,seq0000015:0.13713)0.609:0.00215,seq0000013:0.02032,(seq0000011:0.00014,seq0000012:0.00014)0.766:0.00015);
"""

exp_tree = """((seq0000014:0.08408,seq0000015:0.13713,7:0.02027):0.00215,seq0000013:0.02032,(seq0000011:0.00014,seq0000012:0.00014,6:0.02027):0.00015):0.0;"""

if __name__ == '__main__':
    main()
