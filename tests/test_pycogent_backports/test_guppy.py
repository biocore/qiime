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
from qiime.pycogent_backports.guppy import Guppy, \
                                        build_tree_from_json_using_params
from cogent.app.util import ApplicationError,get_tmp_filename
from cogent.core.tree import PhyloNode
from cogent.parse.tree import DndParser

class Genericguppy(TestCase):

    def setUp(self):
        '''setup the files for testing guppy'''
        
        # create a list of files to cleanup
        self._paths_to_clean_up = []
        self._dirs_to_clean_up = []
        
        # get a tmp filename to use
        basename=splitext(get_tmp_filename())[0]
        
        # create and write out RAxML stats file
        self.json_fname=basename+'.json'
        json_out=open(self.json_fname,'w')
        json_out.write(JSON_RESULT)
        json_out.close()
        self._paths_to_clean_up.append(self.json_fname)

    def tearDown(self): 
        """cleans up all files initially created"""
        # remove the tempdir and contents
        map(remove,self._paths_to_clean_up)
        map(rmdir,self._dirs_to_clean_up)

class guppyTests(Genericguppy):
    """Tests for the guppy application controller"""

    def test_guppy(self):
        """Base command-calls"""
        
        app=Guppy()
        
        self.assertEqual(app.BaseCommand, \
                         ''.join(['cd "',getcwd(),'/"; ','guppy']))
        
        app.Parameters['--help'].on()
        self.assertEqual(app.BaseCommand, \
                         ''.join(['cd "',getcwd(),'/"; ','guppy --help']))
    
    def test_change_working_dir(self):
        """Change working dir"""
        
        # define working directory for output
        working_dir='/tmp/Guppy'
        self._dirs_to_clean_up.append(working_dir)
        
        app = Guppy(WorkingDir=working_dir)
        
        self.assertEqual(app.BaseCommand, \
                       ''.join(['cd "','/tmp/Guppy','/"; ','guppy']))
                    
    def test_build_tree_from_alignment_using_params(self):
        """Builds a tree from a json file"""
        
        # define working directory for output
        outdir='/tmp/'
        
        # set params
        params={}
        
        params["tog"] = None

        # build tree
        tree = build_tree_from_json_using_params(self.json_fname,
                                                      output_dir=outdir,
                                                      params=params)
        
        self.assertEqual(tree.getNewick(),
                    DndParser(TREE_RESULT, constructor=PhyloNode).getNewick())
        
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
    "guppy -t %s -r %s -s %s --out-dir \/tmp %s"
  }, "version": 1, "fields":
  ["edge_num", "likelihood", "like_weight_ratio", "distal_length",
    "pendant_length"
  ]
}
""".replace('\n','').replace(' ','')

TREE_RESULT="""\
((((seq0000004:0.035395,seq0000007:6e-06)\
:0.029109,seq0000006:6e-06)\
:0.019576,seq0000005:0.13713)\
0.609:0.00215,seq0000003:0.02032,(seq0000001:0.00014,seq0000002:0.00014)\
0.766:0.00015)\
:0;
"""

if __name__ == '__main__':
    main()
