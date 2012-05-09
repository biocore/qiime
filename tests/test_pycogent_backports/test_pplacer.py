#!/bin/env python

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"

from os import getcwd, remove, rmdir, mkdir
from os.path import splitext
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from random import randint
from qiime.pycogent_backports.pplacer import Pplacer, insert_sequences_into_tree
from cogent.app.util import ApplicationError,get_tmp_filename
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.tree import PhyloNode
from cogent.core.moltype import RNA,DNA
from StringIO import StringIO
from cogent.core.alignment import Alignment

class Genericpplacer(TestCase):

    def setUp(self):
        '''setup the files for testing pplacer'''
        
        # create a list of files to cleanup
        self._paths_to_clean_up = []
        self._dirs_to_clean_up = []
        
        # get a tmp filename to use
        basename=splitext(get_tmp_filename())[0]
        
        # create and write out RAxML stats file
        self.stats_fname=basename+'.stats'
        stats_out=open(self.stats_fname,'w')
        stats_out.write(RAXML_STATS)
        stats_out.close()
        self._paths_to_clean_up.append(self.stats_fname)
        
        # create and write out reference sequence file
        self.refseq_fname=basename+'_refseqs.fasta'
        refseq_out=open(self.refseq_fname,'w')
        refseq_out.write(REF_SEQS)
        refseq_out.close()
        self._paths_to_clean_up.append(self.refseq_fname)
        
        # create and write out query sequence file
        self.query_fname=basename+'_queryseqs.fasta'
        query_out=open(self.query_fname,'w')
        query_out.write(QUERY_SEQS)
        query_out.close()
        self._paths_to_clean_up.append(self.query_fname)
        
        # create and write out starting tree file
        self.tree_fname=basename+'.tre'
        tree_out=open(self.tree_fname,'w')
        tree_out.write(REF_TREE)
        tree_out.close()
        self._paths_to_clean_up.append(self.tree_fname) 

    def writeTmp(self, outname):
        """Write data to temp file"""
        t = open(outname, "w+")
        t.write(PHYLIP_FILE)
        t.close()

    #
    def tearDown(self): 
        """cleans up all files initially created"""
        # remove the tempdir and contents
        map(remove,self._paths_to_clean_up)
        map(rmdir,self._dirs_to_clean_up)

class pplacerTests(Genericpplacer):
    """Tests for the pplacer application controller"""
    
    def test_pplacer(self):
        """Base command-calls"""
        
        app=Pplacer()
        
        self.assertEqual(app.BaseCommand, \
                         ''.join(['cd "',getcwd(),'/"; ','pplacer']))
        
        app.Parameters['--help'].on()
        self.assertEqual(app.BaseCommand, \
                         ''.join(['cd "',getcwd(),'/"; ','pplacer --help']))
    
    def test_change_working_dir(self):
        """Change working dir"""
        
        working_dir='/tmp/Pplacer'
        self._dirs_to_clean_up.append(working_dir)
        
        # define working directory for output
        app = Pplacer(WorkingDir=working_dir)
        
        self.assertEqual(app.BaseCommand, \
                       ''.join(['cd "','/tmp/Pplacer','/"; ','pplacer']))


    def test_insert_sequences_into_tree(self):
        """Inserts sequences into Tree"""
        
        params={}
        # generate temp filename for output
        params["-r"] = self.refseq_fname
        params["-t"] = self.tree_fname
        params["-s"] = self.stats_fname
        params["--out-dir"] = "/tmp"
        
        aln_ref_query=MinimalFastaParser(StringIO(QUERY_SEQS))
        aln = Alignment(aln_ref_query)
        seqs, align_map = aln.toPhylip()
        tree = insert_sequences_into_tree(seqs, DNA, params=params,
                                          write_log=False)
        
        # rename tips back to query names
        for node in tree.tips():
            if node.Name in align_map:
                node.Name = align_map[node.Name]
        
        self.assertEqual(tree.getNewick(with_distances=True), RESULT_TREE)
        
        
JSON_RESULT="""\
{"tree":
  "((seq0000004:0.08408[0],seq0000005:0.13713[1])0.609:0.00215[2],seq0000003:0.02032[3],(seq0000001:0.00014[4],seq0000002:0.00014[5])0.766:0.00015[6]):0[7];",
  "placements":
  [
    {"p":
      [[0, -113.210938, 0.713818, 0.064504, 0.000006],
        [1, -114.929894, 0.127954, 0.137122, 0.000007],
        [2, -114.932766, 0.127587, 0.000008, 0.000006],
        [6, -117.743534, 0.007675, 0.000141, 0.027211],
        [3, -117.743759, 0.007674, 0.020310, 0.027207],
        [4, -117.747386, 0.007646, 0.000131, 0.027266],
        [5, -117.747396, 0.007646, 0.000131, 0.027266]
      ], "n": ["seq0000006"]
    },
    {"p": [[0, -113.476305, 1.000000, 0.035395, 0.000006]], "n":
      ["seq0000007"]
    }
  ], "metadata":
  {"invocation":
    "pplacer -t %s -r %s -s %s --out-dir \/tmp %s"
  }, "version": 1, "fields":
  ["edge_num", "likelihood", "like_weight_ratio", "distal_length",
    "pendant_length"
  ]
}
""".replace('\n','').replace(' ','')


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

RESULT_TREE="""((((seq0000014:0.0353946,7:6.11352e-06):0.0291093,6:6.11352e-06):0.019576,seq0000015:0.13713)0.609:0.00215,seq0000013:0.02032,(seq0000011:0.00014,seq0000012:0.00014)0.766:0.00015);"""

RAXML_STATS="""


This is RAxML version 7.2.6 released by Alexandros Stamatakis in February 2010.

With greatly appreciated code contributions by:
Andre Aberer (TUM)
Simon Berger (TUM)
John Cazes (TACC)
Michael Ott (TUM)
Nick Pattengale (UNM)
Wayne Pfeiffer (SDSC)


Alignment has 18 distinct alignment patterns

Proportion of gaps and completely undetermined characters in this alignment: 4.80%

RAxML rapid hill-climbing mode

Using 1 distinct models/data partitions with joint branch length optimization


Executing 1 inferences on the original alignment using 1 distinct randomized MP trees

All free model parameters will be estimated by RAxML
ML estimate of 25 per site rate categories

Likelihood of final tree will be evaluated and optimized under GAMMA

GAMMA Model parameters will be estimated up to an accuracy of 0.1000000000 Log Likelihood units

Partition: 0
Alignment Patterns: 18
Name: No Name Provided
DataType: DNA
Substitution Matrix: GTR




RAxML was called as follows:

raxmlHPC -m GTRCAT -s test_raxml.phy -n results 


Inference[0]: Time 0.072128 CAT-based likelihood -85.425107, best rearrangement setting 2
alpha[0]: 1.000000 rates[0] ac ag at cg ct gt: 0.000017 0.037400 0.859448 1.304301 0.000017 1.000000 


Conducting final model optimizations on all 1 trees under GAMMA-based models ....

Inference[0] final GAMMA-based Likelihood: -107.575676 tree written to file /home/RAxML_result.results


Starting final GAMMA-based thorough Optimization on tree 0 likelihood -107.575676 .... 

Final GAMMA-based Score of best tree -107.575676

Program execution info written to /home/RAxML_info.results
Best-scoring ML tree written to: /home/RAxML_bestTree.results

Overall execution time: 0.078965 secs or 0.000022 hours or 0.000001 days
"""

if __name__ == '__main__':
    main()
