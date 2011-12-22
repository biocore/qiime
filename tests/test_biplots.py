#!/usr/bin/env python
# File created on 1 Apr 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

import qiime.biplots as bp
import numpy as np

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import get_random_directory_name
from StringIO import StringIO
from qiime.pycogent_backports.rich_otu_table import nparray_to_ll_mat

class BiplotTests(TestCase):
    
    def setUp(self):
        pass

    def tearDown(self):
        pass
        
    def test_get_taxa(self):
        rand_fname = get_random_directory_name(suppress_mkdir=True)
        rand_fname += '_tmp.txt'
        fout = open(rand_fname,'w')
        lines = ['#Full OTU Counts', \
                     'Taxon\tA\tB\tC', \
                     'Root;Bacteria;Acidobacteria\t0.1\t0.2\t0.3', \
                     'Root;Bacteria;TM7\t0.05\t0.0\t0.3' \
                ]
        fout.write('\n'.join(lines))
        fout.close()
        table = '{"rows": [{"id": "Root;Bacteria;Acidobacteria", "metadata": null}, {"id": "Root;Bacteria;TM7", "metadata": null}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 0.10000000000000001], [0, 1, 0.20000000000000001], [0, 2, 0.29999999999999999], [1, 0, 0.050000000000000003], [1, 2, 0.29999999999999999]], "columns": [{"id": "A", "metadata": null}, {"id": "B", "metadata": null}, {"id": "C", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2579", "matrix_type": "sparse", "shape": [2, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-21T23:50:24.795846", "type": "OTU table", "id": null, "matrix_element_type": "float"}'

        otu_ids = ['Root;Bacteria;Acidobacteria','Root;Bacteria;TM7']
        otu_data = nparray_to_ll_mat(np.array([[0.1,0.3],[0.05,0.3]]))

        res = bp.get_taxa(StringIO(table),sample_ids_kept=['A','C'])
        
        self.assertEqual(res[0],otu_ids)
        self.assertEqual(res[1]._data.items(), otu_data.items())
        
    def test_get_taxa_coords(self):
        otu_table = np.array([  [2,0,0,1],
                                [1,1,1,1],
                                [0,2,2,1]],float)
        sample_names = list('WXYZ')
        otu_names = list('abc')
    
        res = bp.get_taxa_coords(otu_table, [.4,.2,.1,.9])
        otu_coords= range(3)
        otu_coords[0] = .4*2/3 + .9*1/3
        otu_coords[1] = .4*1/4 + .2*1/4 + .1*1/4 + .9*1/4
        otu_coords[2] = .4*0/5 + .2*2/5 + .1*2/5 + .9*1/5
        self.assertFloatEqual(res, otu_coords)
    
    def test_get_taxa_prevalence(self):
        otu_table = np.array([  [2,0,0,1],
                                [1,1,1,1],
                                [0,0,0,0]],float)
        sample_weights = [3,1,1,2]
        res = bp.get_taxa_prevalence(otu_table)
        # print res
        # self.assertFloatEqual(res, np.array([(2/3) + 1/2, 1/3+1+1+1/2, 0])/4) 
        self.assertFloatEqual(res, np.array([(2/3) + 1/2, 1/3+1+1+1/2, 0])/4\
            * 4/(2.5+1/3))                    
        otu_table = np.array([  [2,0,0,1],
                                [1,1,1,1],
                                [0,2,2,1]],float)
        res = bp.get_taxa_prevalence(otu_table)
        # print res
        # self.assertFloatEqual(res, np.array([3,4,5])/12) # if no normalize
        self.assertFloatEqual(res, [0,.5,1])
        
    def test_remove_rare_taxa(self):
        otu_table = np.array([  [2,0,0,1],
                                [1,1,1,1],
                                [0,2,2,1]],float)
        taxdata = {}
        taxdata['prevalence'] = np.array([0,.5,1])
        taxdata['counts'] = otu_table
        taxdata['lineages'] = np.array(['A','B','C'])
        bp.remove_rare_taxa(taxdata,nkeep=2)
        self.assertFloatEqual(taxdata['counts'], otu_table[1:3,:])
        self.assertFloatEqual(taxdata['prevalence'], np.array([.5,1]))
        self.assertFloatEqual(taxdata['lineages'], np.array(['B','C']))

    def test_scale_taxa_data_matrix(self):
        coord = np.array([  [1,4,7,0],
                            [2,5,8,1],
                            [3,6,9,2]],float)
        taxdata = {}
        taxdata['prevalence'] = np.array([0,.5,1])
        taxdata['coord'] = coord
        taxdata['lineages'] = np.array(['Root;A','Root;B','Root;C'])
        pct_var = np.array([100,10,1],dtype="float")

        # with scaling
        res = bp.make_mage_taxa(taxdata,3,pct_var,scaled=True,scalars=None,\
                                    radius=1,
                                    min_taxon_radius=10, max_taxon_radius=20,
                                    taxon_alpha=.7)
        self.assertEqual(res, taxa_mage_scale)

        # without scaling
        res = bp.make_mage_taxa(taxdata,3,pct_var,scaled=False,scalars=None,\
                                    radius=1,
                                    min_taxon_radius=10, max_taxon_radius=20,
                                    taxon_alpha=.7)
        self.assertEqual(res, taxa_mage_no_scale)

    def test_make_biplot_scores_output(self):
        """make_biplot_scores_output correctly formats biplot scores"""
        taxa = {}
        taxa['lineages'] = list('ABC')
        taxa['coord'] = np.array([  [2.1,0.2,0.2,1.4],
                                 [1.1,1.2,1.3,1.5],
                                 [-.3,-2,2.5,1.9]],float)
        res = bp.make_biplot_scores_output(taxa)
        exp = ['#Taxon\tpc0\tpc1\tpc2\tpc3',
               'A\t2.1\t0.2\t0.2\t1.4',
               'B\t1.1\t1.2\t1.3\t1.5',
               'C\t-0.3\t-2.0\t2.5\t1.9',
              ]
        self.assertEqual(res, exp)
    
taxa_mage_no_scale = [\
'@group {Taxa (n=3)} collapsible', \
'@balllist color=white radius=10.0 alpha=0.7 dimension=3 master={taxa_points} nobutton', \
'{A} 1.0 4.0 7.0', \
'@labellist color=white radius=10.0 alpha=0.7 dimension=3 master={taxa_labels} nobutton', \
'{A} 1.0 4.0 7.0', \
'@balllist color=white radius=15.0 alpha=0.7 dimension=3 master={taxa_points} nobutton', \
'{B} 2.0 5.0 8.0', \
'@labellist color=white radius=15.0 alpha=0.7 dimension=3 master={taxa_labels} nobutton', \
'{B} 2.0 5.0 8.0', \
'@balllist color=white radius=20.0 alpha=0.7 dimension=3 master={taxa_points} nobutton', \
'{C} 3.0 6.0 9.0', \
'@labellist color=white radius=20.0 alpha=0.7 dimension=3 master={taxa_labels} nobutton', \
'{C} 3.0 6.0 9.0']

taxa_mage_scale = [\
'@group {Taxa (n=3)} collapsible', \
'@balllist color=white radius=10.0 alpha=0.7 dimension=3 master={taxa_points} nobutton', \
'{A} 1.0 0.4 0.07', \
'@labellist color=white radius=10.0 alpha=0.7 dimension=3 master={taxa_labels} nobutton', \
'{A} 1.0 0.4 0.07', \
'@balllist color=white radius=15.0 alpha=0.7 dimension=3 master={taxa_points} nobutton', \
'{B} 2.0 0.5 0.08', \
'@labellist color=white radius=15.0 alpha=0.7 dimension=3 master={taxa_labels} nobutton', \
'{B} 2.0 0.5 0.08', \
'@balllist color=white radius=20.0 alpha=0.7 dimension=3 master={taxa_points} nobutton', \
'{C} 3.0 0.6 0.09', \
'@labellist color=white radius=20.0 alpha=0.7 dimension=3 master={taxa_labels} nobutton', \
'{C} 3.0 0.6 0.09']

if __name__ == "__main__":
    main()
