#!/usr/bin/env python
#file test_make_3d_plots.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Jesse Stombaugh", "Dan Knights"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Pre-release"

from numpy import array, nan
from StringIO import StringIO
from os.path import exists
from cogent.util.unit_test import TestCase, main
from os import remove
from random import choice, randrange
import shutil
from qiime.make_3d_plots import (make_3d_plots,scale_pc_data_matrix,
                                    auto_radius,make_mage_output,
                                    get_map,get_coord,create_dir,
                                    _make_path,combine_map_label_cols,
                                    process_colorby, linear_gradient,
                                    natsort, make_color_dict,
                                 process_custom_axes, get_custom_coords,
                                 remove_nans,scale_custom_coords,
                                 remove_unmapped_samples)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        self.data={}
        self.data['coord']=[['Sample1','Sample2'],array([[-0.2,0.07],\
                            [-0.04,0.2]]),array([0.7,0.6]),\
                            array([25.00,30.00])]
        self.data['map']=[['#Sample-ID','Day'],['Sample1','Day1'],['Sample2',\
                          'Day1']]
                          
        self.coord_header=["Sample1","Sample2","Sample3"]
        self.coords=array([[-0.219044992,0.079674486,0.09233683],[-0.042258081,\
                       0.000204041,0.024837603],[0.080504323,-0.212014503,\
                       -0.088353435]])
        self.groups={}
        self.groups['Day1']=['Sample1','Sample2','Sample3']
        self.colors={}
        self.colors['Day1']='blue'
        self.pct_var=array([25.00,30.00,35.00])
        self.coord_tups = [("1", "2"), ("3", "2"), ("1", "3")]
        self.colors={"Day1":"blue"}
        self.filename='test_pca.txt'
        self.dir_path='/tmp/' 
        self.prefs={}
        self.prefs['Sample']={}   
        self.prefs['Sample']['column']="Day"
        self.mapping=[["Sample-ID","Day","Type"],["Sample1","Day1","Soil"],\
                      ["Sample2","Day1","Soil"],["Sample3","Day1","Soil"]]
        self.mapping2=[["Sample-ID","Day","Type","Height","Weight"],\
                               ["Sample1","Day1","Soil","10","60"],\
                               ["Sample2","Day1","Soil","20","55"],\
                               ["Sample3","Day1","Soil","30","50"]]
        self.axis_names = 'Height,Weight'
        self._paths_to_clean_up = []
        self._dir_to_clean_up = ''

    def tearDown(self):
        map(remove,self._paths_to_clean_up)
        if self._dir_to_clean_up != '':
            shutil.rmtree(self._dir_to_clean_up)

    def test_natsort(self):
        """natsort should perform numeric comparisons on strings"""
        s = 'sample1 sample2 sample11 sample12'.split()
        self.assertEqual(natsort(s), 
            'sample1 sample2 sample11 sample12'.split())

    def test_make_3d_plots(self):
        """make_3d_plots: main script to create kinemage and html file"""
        obs_kin=make_3d_plots(self.coord_header,self.coords,self.pct_var, \
                          self.mapping,self.prefs)

        self.assertEqual(obs_kin,exp_kin_full)

        # test with custom axes
        custom_axes = ['Height','Weight']
        coord_data = array([[10,60,-0.219044992,0.079674486,0.09233683],
                           [20,55,-0.042258081, 0.000204041,0.024837603],
                           [30,50,0.080504323,-0.212014503,-0.088353435]])
        coords = [self.coord_header, coord_data]
        scale_custom_coords(custom_axes,coords) 
        obs_kin=make_3d_plots(self.coord_header,coords[1],self.pct_var, \
                          self.mapping2,self.prefs, custom_axes)
        self.assertEqual(obs_kin, exp_kin_full_axes)
    
    def test_scale_pc_data_matrix(self):
        """scale_pc_data_matrix: Scales the pc data for use in the 3d plots"""
        exp=array([[-1.56460709e-01,6.82924166e-02,9.23368300e-02],\
                   [-3.01843436e-02,1.74892286e-04,2.48376030e-02],\
                   [5.75030879e-02,-1.81726717e-01,-8.83534350e-02]])

        obs=scale_pc_data_matrix(self.coords, self.pct_var)

        self.assertFloatEqual(obs,exp)
  
    def test_auto_radius(self):
        """auto_radius: determines the radius for the circles in the plot"""
        exp=array([0.00299549315])

        obs=auto_radius(self.coords)

        self.assertFloatEqual(obs,exp)
    
    def test_make_mage_output(self):
        """make_mage_output: Create kinemage string given the data"""
        # test without custom axes
        obs_kin=make_mage_output(self.groups,self.colors,self.coord_header,\
                                 self.coords,self.pct_var)
        self.assertEqual(obs_kin,exp_kin_partial)

        # test with custom axes
        custom_axes = ['Height','Weight']
        coord_data = array([[10,60,-0.219044992,0.079674486,0.09233683],
                           [20,55,-0.042258081, 0.000204041,0.024837603],
                           [30,50,0.080504323,-0.212014503,-0.088353435]])
        coords = [self.coord_header, coord_data]
        scale_custom_coords(custom_axes,coords)
        obs_kin=make_mage_output(self.groups,self.colors,self.coord_header,\
                                 coords[1],self.pct_var, custom_axes)
        self.assertEqual(obs_kin, exp_kin_partial_axes)

    def test_combine_map_label_cols(self):
        """combine_map_label_cols: Combine two or more columns from the \
mapping file"""
        self.combinecolorby=['Day','Type']

        exp=[["Sample-ID","Day","Type","DayType"],\
             ["Sample1","Day1","Soil","Day1Soil"],\
             ["Sample2","Day1","Soil","Day1Soil"],\
             ["Sample3","Day1","Soil","Day1Soil"]]
        obs=combine_map_label_cols(self.combinecolorby,self.mapping)

        self.assertEqual(obs,exp)
    
    def test_create_dir(self):
        """create_dir: creates a directory where the kinemage is stored"""

        alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
        alphabet += alphabet.lower()
        alphabet += "01234567890"

        random_dir_name=''.join([choice(alphabet) for i in range(10)])
        foldername = '/tmp/'+random_dir_name+'/'
            
        self._dir_to_clean_up = foldername
        
        obs=create_dir(foldername,'')
        
        self.assertEqual(obs,foldername)
        self.assertTrue(exists(foldername),'The file was not created in \
the appropriate location')
        
    def test_process_colorby(self):
        """process_colorby: parses the cmd line and determines which columns \
from mapping file to color by"""
        self.colorby='Day'
        exp1={}
        exp1['0']={'column':'Day'}
        obs1,obs2=process_colorby(self.colorby,self.data)

        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,self.data)
    
    def test_process_custom_axes(self):
        """process_custom_axes: Parses the custom_axes \
option from the command line"""
        exp = ['Height','Weight']
        obs=process_custom_axes(self.axis_names)
        self.assertEqual(obs,exp)

    def test_get_custom_coords(self):
        """get_custom_coords: Gets custom axis coords from the mapping file."""
        exp = 1
        custom_axes = ['Height','Weight']
        coords = [self.coord_header, self.coords]
        get_custom_coords(custom_axes, self.mapping2, coords)
        exp = array([[10,60,-0.219044992,0.079674486,0.09233683],
                           [20,55,-0.042258081, 0.000204041,0.024837603],
                           [30,50,0.080504323,-0.212014503,-0.088353435]])
        self.assertEqual(coords[1],exp)

    def test_scale_custom_coords(self):
        """scale_custom_coords: \
Scales custom coordinates to match min/max of PC1"""        
        custom_axes = ['Height','Weight']
        coord_data = array([[10,60,-0.219044992,0.079674486,0.09233683],
                           [20,55,-0.042258081, 0.000204041,0.024837603],
                           [30,50,0.080504323,-0.212014503,-0.088353435]])
        coords = [self.coord_header, coord_data]
        scale_custom_coords(custom_axes,coords)
        # calculate results
        mn = coord_data[2,].min()
        mx = coord_data[2,].max()
        h = array([10.0,20.0,30.0])
        h = (h-min(h))/(max(h)-min(h))
        h = h * (mx-mn) + mn
        w = array([60.0,55.0,50.0])
        w = (w-min(w))/(max(w)-min(w))
        w = w * (mx-mn) + mn
        exp = array([[h[0],w[0],-0.219044992,0.079674486,0.09233683],
                           [h[1],w[1],-0.042258081, 0.000204041,0.024837603],
                           [h[2],w[2],0.080504323,-0.212014503,-0.088353435]])
        self.assertEqual(coords[1],exp)

    def test_remove_nans(self):
        """remove_nans: Deletes any samples with NANs in their coordinates"""
        coord_data = array([[10,60,-0.219044992,0.079674486,0.09233683],
                           [20,55,-0.042258081, nan,0.024837603],
                           [30,50,0.080504323,-0.212014503,-0.088353435]])
        coords = [self.coord_header, coord_data]
        remove_nans(coords)

        exp_header = ["Sample1","Sample3"]
        exp_coords = array([[10,60,-0.219044992,0.079674486,0.09233683],
                           [30,50,0.080504323,-0.212014503,-0.088353435]])
        self.assertEqual(coords[0],exp_header)
        self.assertEqual(coords[1],exp_coords)


    def test_remove_unmapped_samples(self):
        """remove_unmapped_samples: \
Removes any samples not present in mapping file"""
        coord_data = array([[10,60,-0.219044992,0.079674486,0.09233683],
                           [20,55,-0.042258081, nan,0.024837603],
                           [30,50,0.080504323,-0.212014503,-0.088353435]])
        coords = [self.coord_header, coord_data]
        # mapping without sample2
        mapping=[["Sample-ID","Day","Type","Height","Weight"],\
                               ["Sample1","Day1","Soil","10","60"],\
                               ["Sample3","Day1","Soil","30","50"]]
        remove_unmapped_samples(mapping, coords)
        exp_header = ["Sample1","Sample3"]
        exp_coords = array([[10,60,-0.219044992,0.079674486,0.09233683],
                           [30,50,0.080504323,-0.212014503,-0.088353435]])
        self.assertEqual(coords[0],exp_header)
        self.assertEqual(coords[1],exp_coords)


    def test__make_path(self):
        """_make_path: Combine directory and filename into string"""
        obs=_make_path([self.dir_path,self.filename])
        exp='/tmp/test_pca.txt/'

        self.assertEqual(obs,exp)

    def test_make_linear_gradient(self):
        """make_linear_gradient: returns linear gradient of colors"""
        self.assertEqual(linear_gradient([0,1,2],[1,0,0],5),
            [[0,1,2],
             [0.25,0.75,1.5],
             [0.5,0.5,1],
             [0.75,0.25,0.5],
             [1,0,0]])

    def test_make_color_dict(self):
        """make_color_dict: returns dict of named colors"""
        self.assertEqual(make_color_dict('red',(0,100,100),'white',(0,0,100),3),
            {   'redtowhite3_0':[0,100,100],
                'redtowhite3_1':[0,50,100],
                'redtowhite3_2':[0,0,100],
            })
        
exp_kin_full=\
['@kinemage {Day_unscaled}', '@dimension {PC1} {PC2} {PC3}', \
'@dimminmax -0.219044992 0.080504323 -0.212014503 0.079674486 -0.088353435 0.09233683', \
'@master {points}', '@master {labels}', '@hsvcolor {aqua} 180.0 100.0 100.0', \
'@hsvcolor {blue} 240.0 100.0 100.0', \
'@hsvcolor {fuchsia} 300.0 100.0 100.0', '@hsvcolor {gray} 300.0 0.0 50.2', \
'@hsvcolor {green} 120.0 100.0 50.2', '@hsvcolor {lime} 120.0 100.0 100.0', \
'@hsvcolor {maroon} 0.0 100.0 50.2', \
'@hsvcolor {olive} 60.0 100.0 50.2', '@hsvcolor {purple} 300.0 100.0 50.2', \
'@hsvcolor {red} 0.0 100.0 100.0', '@hsvcolor {silver} 0.0 0.0 75.3', \
'@hsvcolor {teal} 180.0 100.0 50.2', \
'@hsvcolor {yellow} 60.0 100.0 100.0', '@hsvcolor {white} 180.0 0.0 100.0', \
'@group {Day1 (n=3)} collapsible', \
'@balllist color=blue radius=0.00299549315 alpha=0.75 dimension=3 master={points} nobutton', \
'{Sample1} -0.219044992 0.079674486 0.09233683\n{Sample2} -0.042258081 0.000204041 0.024837603\n{Sample3} 0.080504323 -0.212014503 -0.088353435', \
'@labellist color=blue radius=0.00299549315 alpha=0.75 dimension=3 master={labels} nobutton', \
'{Sample1} -0.219044992 0.079674486 0.09233683\n{Sample2} -0.042258081 0.000204041 0.024837603\n{Sample3} 0.080504323 -0.212014503 -0.088353435', \
'@group {axes} collapsible', '@vectorlist {PC1 line} dimension=3 on', \
'-0.2299972416 -0.22261522815 -0.09277110675 white', \
'0.08452953915 -0.22261522815 -0.09277110675 white', \
'@labellist {PC1 (25%)} dimension=3 on', \
'{PC1 (25%)}0.0887560161075 -0.22261522815 -0.09277110675 white', \
'@vectorlist {PC2 line} dimension=3 on', \
'-0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 0.0836582103 -0.09277110675 white', \
'@labellist {PC2 (30%)} dimension=3 on', \
'{PC2 (30%)}-0.2299972416 0.087841120815 -0.09277110675 white', \
'@vectorlist {PC3 line} dimension=3 on', \
'-0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 -0.22261522815 0.0969536715 white', \
'@labellist {PC3 (35%)} dimension=3 on', \
'{PC3 (35%)}-0.2299972416 -0.22261522815 0.101801355075 white', \
'@kinemage {Day_scaled}', '@dimension {PC1} {PC2} {PC3}', \
'@dimminmax -0.156460708571 0.0575030878571 -0.181726716857 0.0682924165714 -0.088353435 0.09233683', \
'@master {points}', '@master {labels}', '@hsvcolor {aqua} 180.0 100.0 100.0', \
'@hsvcolor {blue} 240.0 100.0 100.0', \
'@hsvcolor {fuchsia} 300.0 100.0 100.0', '@hsvcolor {gray} 300.0 0.0 50.2', \
'@hsvcolor {green} 120.0 100.0 50.2', '@hsvcolor {lime} 120.0 100.0 100.0', \
'@hsvcolor {maroon} 0.0 100.0 50.2', \
'@hsvcolor {olive} 60.0 100.0 50.2', '@hsvcolor {purple} 300.0 100.0 50.2', \
'@hsvcolor {red} 0.0 100.0 100.0', '@hsvcolor {silver} 0.0 0.0 75.3', \
'@hsvcolor {teal} 180.0 100.0 50.2', \
'@hsvcolor {yellow} 60.0 100.0 100.0', '@hsvcolor {white} 180.0 0.0 100.0', \
'@group {Day1 (n=3)} collapsible', \
'@balllist color=blue radius=0.00213963796429 alpha=0.75 dimension=3 master={points} nobutton', \
'{Sample1} -0.156460708571 0.0682924165714 0.09233683\n{Sample2} -0.0301843435714 0.000174892285714 0.024837603\n{Sample3} 0.0575030878571 -0.181726716857 -0.088353435', \
'@labellist color=blue radius=0.00213963796429 alpha=0.75 dimension=3 master={labels} nobutton', \
'{Sample1} -0.156460708571 0.0682924165714 0.09233683\n{Sample2} -0.0301843435714 0.000174892285714 0.024837603\n{Sample3} 0.0575030878571 -0.181726716857 -0.088353435', \
'@group {axes} collapsible', '@vectorlist {PC1 line} dimension=3 on', \
'-0.164283744 -0.1908130527 -0.09277110675 white', \
'0.06037824225 -0.1908130527 -0.09277110675 white', \
'@labellist {PC1 (25%)} dimension=3 on', \
'{PC1 (25%)}0.0633971543625 -0.1908130527 -0.09277110675 white', \
'@vectorlist {PC2 line} dimension=3 on', \
'-0.164283744 -0.1908130527 -0.09277110675 white', \
'-0.164283744 0.0717070374 -0.09277110675 white', \
'@labellist {PC2 (30%)} dimension=3 on', \
'{PC2 (30%)}-0.164283744 0.07529238927 -0.09277110675 white', \
'@vectorlist {PC3 line} dimension=3 on', \
'-0.164283744 -0.1908130527 -0.09277110675 white', \
'-0.164283744 -0.1908130527 0.0969536715 white', '@labellist {PC3 (35%)} dimension=3 on', \
'{PC3 (35%)}-0.164283744 -0.1908130527 0.101801355075 white']

exp_kin_full_axes =\
['@kinemage {Day_unscaled}', '@dimension {Height} {Weight} {PC1} {PC2} {PC3}', \
'@dimminmax -0.219044992 0.080504323 -0.219044992 0.080504323 -0.219044992 0.080504323 -0.212014503 0.079674486 -0.088353435 0.09233683', \
'@master {points}', '@master {labels}', '@hsvcolor {aqua} 180.0 100.0 100.0', '@hsvcolor {blue} 240.0 100.0 100.0', '@hsvcolor {fuchsia} 300.0 100.0 100.0', '@hsvcolor {gray} 300.0 0.0 50.2', '@hsvcolor {green} 120.0 100.0 50.2', '@hsvcolor {lime} 120.0 100.0 100.0', '@hsvcolor {maroon} 0.0 100.0 50.2', '@hsvcolor {olive} 60.0 100.0 50.2', '@hsvcolor {purple} 300.0 100.0 50.2', '@hsvcolor {red} 0.0 100.0 100.0', '@hsvcolor {silver} 0.0 0.0 75.3', '@hsvcolor {teal} 180.0 100.0 50.2', \
'@hsvcolor {yellow} 60.0 100.0 100.0', '@hsvcolor {white} 180.0 0.0 100.0',\
'@group {Day1 (n=3)} collapsible', \
'@balllist color=blue radius=0.00299549315 alpha=0.75 dimension=5 master={points} nobutton', \
'{Sample1} -0.219044992 0.080504323 -0.219044992 0.079674486 0.09233683\n{Sample2} -0.0692703345 -0.0692703345 -0.042258081 0.000204041 0.024837603\n{Sample3} 0.080504323 -0.219044992 0.080504323 -0.212014503 -0.088353435', \
'@labellist color=blue radius=0.00299549315 alpha=0.75 dimension=5 master={labels} nobutton', \
'{Sample1} -0.219044992 0.080504323 -0.219044992 0.079674486 0.09233683\n{Sample2} -0.0692703345 -0.0692703345 -0.042258081 0.000204041 0.024837603\n{Sample3} 0.080504323 -0.219044992 0.080504323 -0.212014503 -0.088353435', \
'@group {axes} collapsible', '@vectorlist {Height line} dimension=5 on', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'0.08452953915 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'@labellist {Height} dimension=5 on', \
'{Height}0.0887560161075 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'@vectorlist {Weight line} dimension=5 on', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 0.08452953915 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'@labellist {Weight} dimension=5 on', \
'{Weight}-0.2299972416 0.0887560161075 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'@vectorlist {PC1 line} dimension=5 on', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 -0.2299972416 0.08452953915 -0.22261522815 -0.09277110675 white', \
'@labellist {PC1 (25%)} dimension=5 on', \
'{PC1 (25%)}-0.2299972416 -0.2299972416 0.0887560161075 -0.22261522815 -0.09277110675 white', \
'@vectorlist {PC2 line} dimension=5 off', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 -0.2299972416 -0.2299972416 0.0836582103 -0.09277110675 white', \
'@labellist {PC2 (30%)} dimension=5 off', \
'{PC2 (30%)}-0.2299972416 -0.2299972416 -0.2299972416 0.087841120815 -0.09277110675 white', \
'@vectorlist {PC3 line} dimension=5 off', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 0.0969536715 white', \
'@labellist {PC3 (35%)} dimension=5 off', \
'{PC3 (35%)}-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 0.101801355075 white', \
'@kinemage {Day_scaled}', '@dimension {Height} {Weight} {PC1} {PC2} {PC3}', \
'@dimminmax -0.156460708571 0.0575030878571 -0.156460708571 0.0575030878571 -0.156460708571 0.0575030878571 -0.181726716857 0.0682924165714 -0.088353435 0.09233683', \
'@master {points}', '@master {labels}', \
'@hsvcolor {aqua} 180.0 100.0 100.0', '@hsvcolor {blue} 240.0 100.0 100.0', '@hsvcolor {fuchsia} 300.0 100.0 100.0', '@hsvcolor {gray} 300.0 0.0 50.2', '@hsvcolor {green} 120.0 100.0 50.2', '@hsvcolor {lime} 120.0 100.0 100.0', '@hsvcolor {maroon} 0.0 100.0 50.2', '@hsvcolor {olive} 60.0 100.0 50.2', '@hsvcolor {purple} 300.0 100.0 50.2', '@hsvcolor {red} 0.0 100.0 100.0', '@hsvcolor {silver} 0.0 0.0 75.3', '@hsvcolor {teal} 180.0 100.0 50.2', '@hsvcolor {yellow} 60.0 100.0 100.0', '@hsvcolor {white} 180.0 0.0 100.0', '@group {Day1 (n=3)} collapsible', \
'@balllist color=blue radius=0.00213963796429 alpha=0.75 dimension=5 master={points} nobutton', \
'{Sample1} -0.156460708571 0.0575030878571 -0.156460708571 0.0682924165714 0.09233683\n{Sample2} -0.0494788103571 -0.0494788103571 -0.0301843435714 0.000174892285714 0.024837603\n{Sample3} 0.0575030878571 -0.156460708571 0.0575030878571 -0.181726716857 -0.088353435', \
'@labellist color=blue radius=0.00213963796429 alpha=0.75 dimension=5 master={labels} nobutton', \
'{Sample1} -0.156460708571 0.0575030878571 -0.156460708571 0.0682924165714 0.09233683\n{Sample2} -0.0494788103571 -0.0494788103571 -0.0301843435714 0.000174892285714 0.024837603\n{Sample3} 0.0575030878571 -0.156460708571 0.0575030878571 -0.181726716857 -0.088353435', \
'@group {axes} collapsible', '@vectorlist {Height line} dimension=5 on', \
'-0.164283744 -0.164283744 -0.164283744 -0.1908130527 -0.09277110675 white', \
'0.06037824225 -0.164283744 -0.164283744 -0.1908130527 -0.09277110675 white', \
'@labellist {Height} dimension=5 on', \
'{Height}0.0633971543625 -0.164283744 -0.164283744 -0.1908130527 -0.09277110675 white', \
'@vectorlist {Weight line} dimension=5 on', \
'-0.164283744 -0.164283744 -0.164283744 -0.1908130527 -0.09277110675 white', \
'-0.164283744 0.06037824225 -0.164283744 -0.1908130527 -0.09277110675 white', \
'@labellist {Weight} dimension=5 on', \
'{Weight}-0.164283744 0.0633971543625 -0.164283744 -0.1908130527 -0.09277110675 white', \
'@vectorlist {PC1 line} dimension=5 on', \
'-0.164283744 -0.164283744 -0.164283744 -0.1908130527 -0.09277110675 white', \
'-0.164283744 -0.164283744 0.06037824225 -0.1908130527 -0.09277110675 white', \
'@labellist {PC1 (25%)} dimension=5 on', \
'{PC1 (25%)}-0.164283744 -0.164283744 0.0633971543625 -0.1908130527 -0.09277110675 white', \
'@vectorlist {PC2 line} dimension=5 off', \
'-0.164283744 -0.164283744 -0.164283744 -0.1908130527 -0.09277110675 white', \
'-0.164283744 -0.164283744 -0.164283744 0.0717070374 -0.09277110675 white', \
'@labellist {PC2 (30%)} dimension=5 off', \
'{PC2 (30%)}-0.164283744 -0.164283744 -0.164283744 0.07529238927 -0.09277110675 white', \
'@vectorlist {PC3 line} dimension=5 off', \
'-0.164283744 -0.164283744 -0.164283744 -0.1908130527 -0.09277110675 white', \
'-0.164283744 -0.164283744 -0.164283744 -0.1908130527 0.0969536715 white', \
'@labellist {PC3 (35%)} dimension=5 off', \
'{PC3 (35%)}-0.164283744 -0.164283744 -0.164283744 -0.1908130527 0.101801355075 white']

exp_kin_partial=\
['@kinemage {_unscaled}', '@dimension {PC1} {PC2} {PC3}', \
'@dimminmax -0.219044992 0.080504323 -0.212014503 0.079674486 -0.088353435 0.09233683', \
'@master {points}', '@master {labels}', '@hsvcolor {aqua} 180.0 100.0 100.0', \
'@hsvcolor {blue} 240.0 100.0 100.0', \
'@hsvcolor {fuchsia} 300.0 100.0 100.0', '@hsvcolor {gray} 300.0 0.0 50.2', \
'@hsvcolor {green} 120.0 100.0 50.2', '@hsvcolor {lime} 120.0 100.0 100.0', \
'@hsvcolor {maroon} 0.0 100.0 50.2', \
'@hsvcolor {olive} 60.0 100.0 50.2', '@hsvcolor {purple} 300.0 100.0 50.2', \
'@hsvcolor {red} 0.0 100.0 100.0', '@hsvcolor {silver} 0.0 0.0 75.3', \
'@hsvcolor {teal} 180.0 100.0 50.2', \
'@hsvcolor {yellow} 60.0 100.0 100.0', '@hsvcolor {white} 180.0 0.0 100.0',\
'@group {Day1 (n=3)} collapsible', \
'@balllist color=blue radius=0.00299549315 alpha=0.75 dimension=3 master={points} nobutton', \
'{Sample1} -0.219044992 0.079674486 0.09233683\n{Sample2} -0.042258081 0.000204041 0.024837603\n{Sample3} 0.080504323 -0.212014503 -0.088353435', \
'@labellist color=blue radius=0.00299549315 alpha=0.75 dimension=3 master={labels} nobutton', \
'{Sample1} -0.219044992 0.079674486 0.09233683\n{Sample2} -0.042258081 0.000204041 0.024837603\n{Sample3} 0.080504323 -0.212014503 -0.088353435', \
'@group {axes} collapsible', '@vectorlist {PC1 line} dimension=3 on', \
'-0.2299972416 -0.22261522815 -0.09277110675 white', \
'0.08452953915 -0.22261522815 -0.09277110675 white', \
'@labellist {PC1 (25%)} dimension=3 on', '{PC1 (25%)}0.0887560161075 -0.22261522815 -0.09277110675 white', '@vectorlist {PC2 line} dimension=3 on', \
'-0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 0.0836582103 -0.09277110675 white', \
'@labellist {PC2 (30%)} dimension=3 on', \
'{PC2 (30%)}-0.2299972416 0.087841120815 -0.09277110675 white', \
'@vectorlist {PC3 line} dimension=3 on', \
'-0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 -0.22261522815 0.0969536715 white', \
'@labellist {PC3 (35%)} dimension=3 on', \
'{PC3 (35%)}-0.2299972416 -0.22261522815 0.101801355075 white']


exp_kin_partial_axes=\
['@kinemage {_unscaled}', '@dimension {Height} {Weight} {PC1} {PC2} {PC3}',\
'@dimminmax -0.219044992 0.080504323 -0.219044992 0.080504323 -0.219044992 0.080504323 -0.212014503 0.079674486 -0.088353435 0.09233683',
'@master {points}', '@master {labels}', '@hsvcolor {aqua} 180.0 100.0 100.0', \
'@hsvcolor {blue} 240.0 100.0 100.0', \
'@hsvcolor {fuchsia} 300.0 100.0 100.0', '@hsvcolor {gray} 300.0 0.0 50.2', \
'@hsvcolor {green} 120.0 100.0 50.2', '@hsvcolor {lime} 120.0 100.0 100.0', \
'@hsvcolor {maroon} 0.0 100.0 50.2', '@hsvcolor {olive} 60.0 100.0 50.2', \
'@hsvcolor {purple} 300.0 100.0 50.2', '@hsvcolor {red} 0.0 100.0 100.0', \
'@hsvcolor {silver} 0.0 0.0 75.3', '@hsvcolor {teal} 180.0 100.0 50.2', \
'@hsvcolor {yellow} 60.0 100.0 100.0', '@hsvcolor {white} 180.0 0.0 100.0',\
'@group {Day1 (n=3)} collapsible',
'@balllist color=blue radius=0.00299549315 alpha=0.75 dimension=5 master={points} nobutton', \
'{Sample1} -0.219044992 0.080504323 -0.219044992 0.079674486 0.09233683\n{Sample2} -0.0692703345 -0.0692703345 -0.042258081 0.000204041 0.024837603\n{Sample3} 0.080504323 -0.219044992 0.080504323 -0.212014503 -0.088353435', \
'@labellist color=blue radius=0.00299549315 alpha=0.75 dimension=5 master={labels} nobutton', \
'{Sample1} -0.219044992 0.080504323 -0.219044992 0.079674486 0.09233683\n{Sample2} -0.0692703345 -0.0692703345 -0.042258081 0.000204041 0.024837603\n{Sample3} 0.080504323 -0.219044992 0.080504323 -0.212014503 -0.088353435', \
'@group {axes} collapsible', '@vectorlist {Height line} dimension=5 on', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'0.08452953915 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'@labellist {Height} dimension=5 on', \
'{Height}0.0887560161075 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'@vectorlist {Weight line} dimension=5 on', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 0.08452953915 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'@labellist {Weight} dimension=5 on', \
'{Weight}-0.2299972416 0.0887560161075 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'@vectorlist {PC1 line} dimension=5 on', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 -0.2299972416 0.08452953915 -0.22261522815 -0.09277110675 white', \
'@labellist {PC1 (25%)} dimension=5 on', \
'{PC1 (25%)}-0.2299972416 -0.2299972416 0.0887560161075 -0.22261522815 -0.09277110675 white', \
'@vectorlist {PC2 line} dimension=5 off', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 -0.2299972416 -0.2299972416 0.0836582103 -0.09277110675 white', \
'@labellist {PC2 (30%)} dimension=5 off', \
'{PC2 (30%)}-0.2299972416 -0.2299972416 -0.2299972416 0.087841120815 -0.09277110675 white', \
'@vectorlist {PC3 line} dimension=5 off', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 -0.09277110675 white', \
'-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 0.0969536715 white', \
'@labellist {PC3 (35%)} dimension=5 off', \
'{PC3 (35%)}-0.2299972416 -0.2299972416 -0.2299972416 -0.22261522815 0.101801355075 white']





#run tests if called from command line
if __name__ == "__main__":
    main()
