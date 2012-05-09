#!/bin/env python

from os import getcwd, remove, rmdir, mkdir
from os.path import splitext
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from random import randint
from qiime.pycogent_backports.raxml_v730 import (Raxml,raxml_alignment,\
                                   build_tree_from_alignment,\
                                   insert_sequences_into_tree)
from cogent.app.util import ApplicationError,get_tmp_filename
from cogent.parse.phylip import get_align_for_phylip
from cogent.core.tree import PhyloNode
from cogent.core.moltype import RNA,DNA
from StringIO import StringIO
from cogent.util.misc import app_path
from subprocess import Popen, PIPE, STDOUT
from cogent.core.alignment import Alignment
import re
from random import choice, randint

__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Micah Hamady", "Rob Knight", "Daniel McDonald","Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Micah Hamady"
__email__ = "Micah Hamady"
__status__ = "Release"

class GenericRaxml(TestCase):

    def setUp(self):
        """Check if Raxml version is supported for this test"""
        acceptable_version = (7,3,0)
        self.assertTrue(app_path('raxmlHPC'),
         "raxmlHPC not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        command = "raxmlHPC -v | grep version"
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[4].strip()
        try:
            version = tuple(map(int,version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported raxmlHPC version. %s is required, but running %s." \
         % ('.'.join(map(str,acceptable_version)), version_string))
        
        
        """Setup data for raxml tests"""
        self.seqs1 = ['ACUGCUAGCUAGUAGCGUACGUA','GCUACGUAGCUAC',
            'GCGGCUAUUAGAUCGUA']
        self.labels1 = ['>1','>2','>3']
        self.lines1 = flatten(zip(self.labels1,self.seqs1))

        self.test_model = "GTRCAT"

        self.align1 = get_align_for_phylip(StringIO(PHYLIP_FILE))

        self.test_fn1 = "/tmp/raxml_test1.txt"
        self.test_fn2 = "raxml_test1.txt"
        self.test_fn1_space = "/tmp/raxml test1.txt"

    def writeTmp(self, outname):
        """Write data to temp file"""
        t = open(outname, "w+")
        t.write(PHYLIP_FILE)
        t.close()


class RaxmlTests(GenericRaxml):
    """Tests for the Raxml application controller"""

    def test_raxml(self):
        """raxml BaseCommand should return the correct BaseCommand"""
        r = Raxml()
        self.assertEqual(r.BaseCommand, \
            ''.join(['cd \"',getcwd(),'/\"; ','raxmlHPC -f d -# 1']))
        r.Parameters['-s'].on('seq.nexus')
        self.assertEqual(r.BaseCommand,\
            ''.join(['cd \"',getcwd(),'/\"; ',\
            'raxmlHPC -f d -s seq.nexus -# 1']))


    def test_raxml_params(self):
        """raxml should raise exception if missing required params"""

        r = Raxml(WorkingDir="/tmp")

        r.SuppressStdout = True
        r.SuppressStderr = True
        # raise error by default
        self.assertRaises(ValueError, r)

        # specify output name 
        r.Parameters['-n'].on("test_name")
        r.Parameters["-p"].on(randint(1,100000))
        self.assertRaises(ApplicationError, r)

        # specify model 
        r.Parameters['-m'].on("GTRCAT")
        self.assertRaises(ApplicationError, r)

        r.Parameters['-s'].on(self.test_fn1)
        self.assertRaises(ApplicationError, r)


        self.writeTmp(self.test_fn1)

        o = r()
        o.cleanUp()

        remove(self.test_fn1)

    
    def test_raxml_from_file(self):
        """raxml should run correctly using filename"""
        r = Raxml(WorkingDir="/tmp")

        r.Parameters['-s'].on(self.test_fn1)
        r.Parameters['-m'].on("GTRCAT")
        r.Parameters['-n'].on("test_me")
        r.Parameters["-p"].on(randint(1,100000))
       
        # test with abs filename
        cur_out = self.test_fn1
        self.writeTmp(cur_out)
        out = r()
        out.cleanUp()
        remove(cur_out)

        # test with rel + working dir 
        r.Parameters['-s'].on(self.test_fn2)
        r.Parameters['-n'].on("test_me2")
        r.Parameters['-w'].on("/tmp/")
        r.Parameters["-p"].on(randint(1,100000))
        self.writeTmp(self.test_fn1)
        out = r()
        out.cleanUp()
        remove(self.test_fn1)

        r.Parameters['-s'].on("\"%s\"" % self.test_fn1_space)
        r.Parameters['-n'].on("test_me3")
        r.Parameters['-w'].on("/tmp/")
        r.Parameters["-p"].on(randint(1,100000))
        #print r.BaseCommand
        self.writeTmp(self.test_fn1_space)
        out = r()
        out.cleanUp()
        remove(self.test_fn1_space)

    def test_raxml_alignment(self):
        """raxml_alignment should work as expected"""
        phy_node, parsimony_phy_node, log_likelihood, total_exec \
            = raxml_alignment(self.align1)

    def test_build_tree_from_alignment(self):
        """Builds a tree from an alignment"""
        
        tree = build_tree_from_alignment(self.align1, RNA, False)
        
        self.assertTrue(isinstance(tree, PhyloNode))
        self.assertEqual(len(tree.tips()), 7)
        self.assertRaises(NotImplementedError, build_tree_from_alignment, \
                          self.align1, RNA, True)
    
    def test_insert_sequences_into_tree(self):
        """Inserts sequences into Tree using params - test handles tree-insertion"""
        
        # generate temp filename for output
        outfname=splitext(get_tmp_filename('/tmp/'))[0]
        
        # create starting tree
        outtreefname=outfname+'.tre'
        outtree=open(outtreefname,'w')
        outtree.write(REF_TREE)
        outtree.close()
        
        # set params for tree-insertion
        params={}
        params["-w"]="/tmp/"
        params["-n"] = get_tmp_filename().split("/")[-1]
        params["-f"] = 'v'
        #params["-G"] = '0.25'
        params["-t"] = outtreefname
        params["-m"] = 'GTRGAMMA'
        
        aln_ref_query=get_align_for_phylip(StringIO(PHYLIP_FILE_DNA_REF_QUERY))
        aln = Alignment(aln_ref_query)
        seqs, align_map = aln.toPhylip()
        
        tree = insert_sequences_into_tree(seqs, DNA, params=params,
                                          write_log=False)
        
        for node in tree.tips():
            removed_query_str=re.sub('QUERY___','',str(node.Name))
            new_node_name=re.sub('___\d+','',str(removed_query_str))
            if new_node_name in align_map:
                node.Name = align_map[new_node_name]

        self.assertTrue(isinstance(tree, PhyloNode))
        self.assertEqual(tree.getNewick(with_distances=True),RESULT_TREE)
        self.assertEqual(len(tree.tips()), 7)
        self.assertRaises(NotImplementedError, build_tree_from_alignment, \
                         self.align1, RNA, True)
                         
        remove(outtreefname)

PHYLIP_FILE= """ 7 50
Species001   UGCAUGUCAG UAUAGCUUUA GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
Species002   UGCAUGUCAG UAUAGCUUUA GUGAAACUGC GAAUGGCUNN UUAAAUCAGU
Species003   UGCAUGUCAG UAUAGCAUUA GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
Species004   UCCAUGUCAG UAUAACUUUG GUGAAACUGC GAAUGGCUCA UUAAAUCAGG
Species005   NNNNNNNNNN UAUAUCUUAU GUGAAACUUC GAAUGCCUCA UUAAAUCAGU
Species006   UGCAUGUCAG UAUAGCUUUG GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
Species007   UGCAUGUCAG UAUAACUUUG GUGAAACUGC GAAUGGCUCA UUAAAUCAGU
""" 


PHYLIP_FILE_DNA_REF_QUERY= """ 7 50
Species001   TGCATGTCAG TATAGCTTTA GTGAAACTGC GAATGGCTCA TTAAATCAGT
Species002   TGCATGTCAG TATAGCTTTA GTGAAACTGC GAATGGCTNN TTAAATCAGT
Species003   TGCATGTCAG TATAGCATTA GTGAAACTGC GAATGGCTCA TTAAATCAGT
Species004   TCCATGTCAG TATAACTTTG GTGAAACTGC GAATGGCTCA TTAAATCAGG
Species005   NNNNNNNNNN TATATCTTAT GTGAAACTTC GAATGCCTCA TTAAATCAGT
Species006   TGCATGTCAG TATAGCTTTG GTGAAACTGC GAATGGCTCA TTAAATCAGT
Species007   TGCATGTCAG TATAACTTTG GTGAAACTGC GAATGGCTCA TTAAATCAGT
"""

REF_TREE="""((seq0000004:0.08408,seq0000005:0.13713)0.609:0.00215,seq0000003:0.02032,(seq0000001:0.00014,seq0000002:0.00014)0.766:0.00015);
"""

RESULT_TREE="""(Species003:0.0194919169324,(Species001:4.34281710439e-07,Species002:4.34281710439e-07):4.34281710439e-07,(((Species006:0.0,Species007:0.0):0.0,Species004:0.0438017433031):0.0438017433031,Species005:0.171345128781):0.00331197405878);"""

if __name__ == '__main__':
    main()
