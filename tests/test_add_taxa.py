#!/usr/bin/env python

"""Tests of code for adding taxa to OTU table"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
#remember to add yourself if you make changes
__credits__ = ["Rob Knight"] 
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
from qiime.add_taxa import (fix_taxonomy_delimiters, 
    rewrite_otu_table_with_taxonomy)
from StringIO import StringIO

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_fix_taxonomy_delimiters(self):
        """fix_taxonomy_delimiters should remove commas from RDP"""
        rdp = { '2483210':['Root,Bacteria','0.940'],
                '2381498':['Root,Bacteria,Firmicutes,"Clostridia",Clostridiales,"Lachnospiraceae"','1.000']}
        self.assertEqual(fix_taxonomy_delimiters(rdp),
            {'2483210':'Root;Bacteria',
             '2381498':'Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae',
             })

    def test_rewrite_otu_table_with_taxonomy(self):
        """rewrite_otu_table_with_taxonomy should add taxonomy string"""
        otu_lines = """#Full OTU Counts
#OTU ID\tPC.354\tPC.355
0\t0\t0
1\t0\t0
2\t0\t0""".splitlines()
        tax_lines = """0 PC.636_424\tRoot;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae\t1
1 PC.481_321\tRoot;Bacteria;Firmicutes;"Clostridia";Clostridiales\t0.89
2 PC.635_886\tRoot;Bacteria\t0.94
""".splitlines()
        result = """#Full OTU Counts
#OTU ID\tPC.354\tPC.355\tConsensus Lineage
0\t0\t0\tRoot;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
1\t0\t0\tRoot;Bacteria;Firmicutes;Clostridia;Clostridiales
2\t0\t0\tRoot;Bacteria
"""
        outfile = StringIO()
        rewrite_otu_table_with_taxonomy(tax_lines, otu_lines, outfile=outfile)
        outfile.seek(0)
        self.assertEqual(outfile.read(), result)

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
