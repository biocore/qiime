#!/usr/bin/env python
# File created on 1 Apr 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

import qiime.biplots as bp
import numpy as np

from os import system
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import get_random_directory_name

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

        otu_ids = ['Root;Bacteria;Acidobacteria','Root;Bacteria;TM7']
        otu_table = np.array([[0.1,0.3],[0.05,0.3]])
        res = bp.get_taxa(rand_fname,sample_ids_kept=['A','C'])
        self.assertEqual(res[0],otu_ids)
        self.assertEqual(res[1],otu_table)

        # remove temporary file
        system('rm %s' %(rand_fname))
        pass
        
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
        coord = np.array([  [1,4,7,0],
                            [2,5,8,1],
                            [3,6,9,2]],float)
        taxdata = {}
        taxdata['coord'] = coord
        taxdata['prevalence'] = np.array([0,.5,1])
        taxdata['counts'] = otu_table
        taxdata['lineages'] = np.array(['A','B','C'])
        bp.remove_rare_taxa(taxdata,nkeep=2)
        
        exp_otu_table = np.array([[0,2,2,1],
                                  [1,1,1,1]],float)
        exp_prevalence = np.array([1,.5])
        exp_lineages = np.array(['C','B'])
        self.assertFloatEqual(taxdata['counts'],     exp_otu_table)
        self.assertFloatEqual(taxdata['prevalence'], exp_prevalence)
        self.assertFloatEqual(taxdata['lineages'],   exp_lineages)

    def test_remove_rare_taxa_not_valid_number(self):
        otu_table = np.array([  [2,0,0,1],
                                [1,1,1,1],
                                [0,2,2,1]],float)
        coord = np.array([  [1,4,7,0],
                            [2,5,8,1],
                            [3,6,9,2]],float)
        taxdata = {}
        taxdata['coord'] = coord
        taxdata['prevalence'] = np.array([0,.5,1])
        taxdata['counts'] = otu_table
        taxdata['lineages'] = np.array(['A','B','C'])
        bp.remove_rare_taxa(taxdata,nkeep=-2)
        
        exp_otu_table = np.array([[0,2,2,1],
                                  [1,1,1,1],
                                  [2,0,0,1]],float)
        exp_prevalence = np.array([1,.5,0])
        exp_lineages = np.array(['C','B','A'])
        self.assertFloatEqual(taxdata['counts'],     exp_otu_table)
        self.assertFloatEqual(taxdata['prevalence'], exp_prevalence)
        self.assertFloatEqual(taxdata['lineages'],   exp_lineages)

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
