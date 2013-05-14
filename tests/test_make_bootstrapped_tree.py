#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["justin kuczynski"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

from qiime.make_bootstrapped_tree import write_pdf_bootstrap_tree
from cogent.util.unit_test import TestCase, main
from cogent.core.tree import PhyloNode
from qiime.parse import parse_newick
from qiime.util import get_tmp_filename
import os


def remove_files(list_of_filepaths,error_on_missing=True):
    missing = []
    for fp in list_of_filepaths:
        try:
            os.remove(fp)
        except OSError:
            missing.append(fp)

    if error_on_missing and missing:
        raise OSError, "Some filepaths were not accessible: %s"\
            % '\t'.join(missing)

class FunctionTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        self._paths_to_clean_up = []
        
    def test_write_pdf_bootstrap_tree(self):
        """ write_pdf_bootstrap_tree should throw no errors"""
        
        tree = parse_newick(
            "((tax7:0.1,tax3:0.2)node0:.98,tax8:.3, tax4:.3)node1:.4",
            PhyloNode)
        bootstraps = {'node0':.7,'node1':.4}
        f = get_tmp_filename(\
         prefix='make_bootstrapped_tree_test',\
         suffix='.pdf',\
         result_constructor=str)
        self._paths_to_clean_up.append(f)
        write_pdf_bootstrap_tree(tree, f, bootstraps)
        assert(os.path.exists(f))
        
    def test_write_pdf_bootstrap_tree_escaped_names(self):
        """ write_pdf_bootstrap_tree functions when newick names are escaped
        
            This test essentially is only checking that no failures arise from
            having escaped strings as nodes in the newick file. Visual inspection
            of the resulting PDFs shows that the coloring is occuring as expected
            but unfortunately there is not a great way to test for this.
        
        """
        
        tree = parse_newick(
            "((tax7:0.1,'tax3':0.2)'no__``!!:o de0':.98,'ta___x8':.3, tax4:.3)node1:.4",
            PhyloNode)
        
        bootstraps = {"no__``!!:o de0":.7,'node1':.4}
        f = get_tmp_filename(\
         prefix='make_bootstrapped_tree_test',\
         suffix='.pdf',\
         result_constructor=str)
        self._paths_to_clean_up.append(f)
        write_pdf_bootstrap_tree(tree, f, bootstraps)
        assert(os.path.exists(f))
        

    def tearDown(self):
        remove_files(self._paths_to_clean_up)
        

#run tests if called from command line
if __name__ == '__main__':
    main()
