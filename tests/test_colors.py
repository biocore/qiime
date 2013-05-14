#!/usr/bin/env python
#file test_colors.py

__author__ = "Rob Knight and Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"

from numpy import array
from StringIO import StringIO
from os.path import exists
from collections import defaultdict
from cogent.util.unit_test import TestCase, main
from os import remove
from random import choice, randrange
import shutil
from qiime.colors import (Color, rgb_tuple_to_hsv, mage_hsv_tuple_to_rgb,
    combine_map_label_cols, process_colorby, 
    linear_gradient, natsort, make_color_dict, color_dict_to_objects,
    iter_color_groups, get_group_colors,
    get_color, color_groups,string_to_rgb,
    get_map,map_from_coords,sample_color_prefs_and_map_data_from_options,
    taxonomy_process_prefs)

class ColorTests(TestCase):
    """Tests of the Color class"""

    def setUp(self):
        """Color init should take name, coords, colorspace"""
        self.black = Color('black', (0,0,0))
        self.red = Color('red', (255,0,0))
        self.pink = Color('pink', (100,0,0))
        self.green = Color('green', (0,255,0))

    def test_str(self):
        """Color init and string representation should give correct result"""
        self.assertEqual(str(self.black), 'black:#000000')
        self.assertEqual(str(self.red), 'red:#ff0000')
        self.assertEqual(str(self.pink), 'pink:#640000')

    def test_toRGB(self):
        """Color toRGB should give correct r, g, b tuple (range 0-255)"""
        self.assertEqual(self.black.toRGB(), (0,0,0))
        self.assertEqual(self.red.toRGB(), (255,0,0))
        self.assertEqual(self.pink.toRGB(), (100,0,0))

    def test_toMage(self):
        """Color toMage should give correct string using h, s, v tuple"""
        self.assertEqual(self.black.toMage(), '@hsvcolor {black} 0.0 0.0 0.0')
        self.assertEqual(self.red.toMage(), '@hsvcolor {red} 0.0 100.0 100.0')
        self.assertEqual(self.pink.toMage(), '@hsvcolor {pink} 0.0 100.0 39.2')
        self.assertEqual(self.green.toMage(), \
            '@hsvcolor {green} 120.0 100.0 100.0')

    def test_toHex(self):
        """Color toHex should give correct hex string"""
        self.assertEqual(self.black.toHex(), '#000000')
        self.assertEqual(self.red.toHex(), '#ff0000')
        self.assertEqual(self.pink.toHex(), '#640000')
        
    def test_toInt(self):
        """Color toHex should give correct hex string"""
        self.assertEqual(self.black.toInt(), 0)
        self.assertEqual(self.red.toInt(), 16711680)
        self.assertEqual(self.pink.toInt(), 6553600)


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
        self.pct_var=array([25.00,30.00,35.00])
        self.coord_tups = [("1", "2"), ("3", "2"), ("1", "3")]
        self.colors={"Day1":"red1"}
        self.filename='test_pca.txt'
        self.dir_path='/tmp/' 
        self.prefs={}
        self.prefs['Sample']={}   
        self.prefs['Sample']['column']="Day"
       
        self.dict=defaultdict(list)
        self.dict['Day1'].append('Sample1')
        self.dict['Day1'].append('Sample2')
        self.dict['Day1'].append('Sample3')
        
        self.labelname=self.prefs['Sample']['column']
        self.mapping=[["Sample-ID","Day","Type"],["Sample1","Day1","Soil"],\
                      ["Sample2","Day1","Soil"],["Sample3","Day1","Soil"]]
        self.data_color_hsv = {
            #'black1':	(0,0,20),
            'red1':	(0,100,100),
            'blue1':	(240,100,100),
            'orange1':	(28,98,95),
            'green1':	(120,100,50.2),
            'purple1':	(302,73,57),
            'yellow1':	(60,100,100),
            'cyan1':	(184, 49, 96),
            'pink1':	(333,37,96),
            'teal1':	(178,42,63),
            'brown1':	(36,89,42),
            'gray1':	(0,0,50.2),
            'lime':	(123,99,96),
            'red2':	(14,51,97),
            'blue2':	(211,42,85),
            'orange2':	(32,46,99),
            'green2':	(142,36,79),
            'purple2':	(269,29,75),
            'yellow2':	(56,40,100),
            #'black2':	(303,100,24),
            'gray2':	(0, 0, 75.3),
            #'teal2':	(192,100,24),
            'red3':	(325,100,93),
            'blue3':	(197,100,100),
            #'purple3':	(271,43,36),
            'brown2':	(33,45,77),
            'green3':	(60,100,50.2),
            'purple4':	(264,75,100),
            #'yellow3':	(60,66,75),
            #'blue4':	(213,45,77),
            'red4':	(348,31,74),
            'teal3':	(180,100,50.2),
            #'brown3':	(60,100,28),
            'red5':	(0,100,50.2),
            'green4':	(81,100,26),
            #'purple5':	(240,100,41),
            'orange3':	(26,100,65)
            #'brown4':	(25,100,20),
            #'red6':	(17,100,63),
            #'purple6':(272,100,44)
        }
        
        self.data_color_order = ['red1', 'blue1', 'orange1', 'green1',\
                    'purple1', 'yellow1', 'cyan1', 'pink1', 'teal1', 'brown1',\
                    'gray1', 'lime', 'red2', 'blue2', 'orange2', 'green2',\
                    'purple2', 'yellow2', 'gray2', 'red3', 'blue3', 'brown2',\
                    'green3', 'purple4', 'red4', 'teal3', 'red5', 'green4',\
                    'orange3']
        
        self._paths_to_clean_up = []
        self._dir_to_clean_up = ''

    def tearDown(self):
        map(remove,self._paths_to_clean_up)
        if self._dir_to_clean_up != '':
            shutil.rmtree(self._dir_to_clean_up)

    def test_string_to_rgb(self):
        """str_to_rgb should accept a hex string and emit tuples on right scale"""
        strgb = string_to_rgb #for convenience
        self.assertEqual(strgb('#000000'), (0,0,0))
        self.assertEqual(strgb('#FFFFFF'), (255,255,255))
        self.assertEqual(strgb('#F0F0F0'), (240,240,240))
        self.assertEqual(strgb('#AFF0AA'), (175, 240, 170))

    def test_rgb_tuple_to_hsv(self):
        """rgb_tuple_to_hsv should accept and emit tuples on right scale"""
        rth = rgb_tuple_to_hsv #for convenience
        self.assertEqual(rth((0,0,0)), (0,0,0))
        self.assertEqual(rth((255,0,0)), (0,100,100))
        self.assertEqual(rth((0,255,0)), (120,100,100))
        self.assertEqual(rth((0,0,255)), (240,100,100))
        self.assertFloatEqual(rth((127,127,127)), (0,0,49.803921568627452))

    def test_mage_hsv_tuple_to_rgb(self):
        """mage_hsv_to_rgb should accept and emit tuples on right scale"""
        htr = mage_hsv_tuple_to_rgb #for convenience
        self.assertEqual(htr((0,0,0)), (0,0,0))
        self.assertEqual(htr((0,100,100)), (255,0,0))
        self.assertEqual(htr((120,100,100)), (0,255,0))
        self.assertEqual(htr((240,100,100)), (0,0,255))
        self.assertFloatEqual(htr((0,0,49.803921568627452)), (127,127,127))

    def test_combine_map_label_cols(self):
        """combine_map_label_cols: Combine two or more columns from the \
mapping file"""
        self.combinecolorby=['Day','Type']
        
        exp=[["Sample-ID","Day","Type","Day&&Type"],\
             ["Sample1","Day1","Soil","Day1Soil"],\
             ["Sample2","Day1","Soil","Day1Soil"],\
             ["Sample3","Day1","Soil","Day1Soil"]]
        
        obs=combine_map_label_cols(self.combinecolorby,self.mapping)

        self.assertEqual(obs,exp)
    
    def test_process_colorby(self):
        """process_colorby: parses the cmd line and determines which columns \
from mapping file to color by"""
        self.colorby='Day'
        exp1={}
        exp1['0']={'column':'Day'}
        obs1,obs2=process_colorby(self.colorby,self.data)

        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,self.data)
    
    def test_linear_gradient(self):
        """linear_gradient: returns linear gradient of colors"""
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

    def test_color_dict_to_objects(self):
        """color_dict_to_objects should return color objects"""
        d = {   'redtowhite3_0':[0,100,100],
                'redtowhite3_1':[0,50,100],
                'redtowhite3_2':[0,0,100],
            }
        res = color_dict_to_objects(d)
        obs = [(str(k), str(v)) for k, v in sorted(res.items())]
        exp = [('redtowhite3_0', 'redtowhite3_0:#ff0000'), 
               ('redtowhite3_1', 'redtowhite3_1:#ff7f7f'), 
               ('redtowhite3_2', 'redtowhite3_2:#ffffff')]
        self.assertEqual(obs, exp)

    def test_iter_color_groups(self):
        """iter_color_groups should iterate over color groups correctly."""

        obs=iter_color_groups(self.mapping,self.prefs)
        obs1=list(obs)
        obs_label=obs1[0][0]
        obs_groups=obs1[0][1]
        obs_colors=obs1[0][2]
        obs_data_colors=obs1[0][3]
        obs_data_color_order=obs1[0][4]
        
        
        data_colors = color_dict_to_objects(self.data_color_hsv)
        
        self.assertEqual(obs_label,self.labelname)
        self.assertEqual(obs_groups,self.dict)
        self.assertEqual(obs_colors,self.colors)
        self.assertEqual(obs_data_colors.keys(),data_colors.keys())
        
        #Need to iterate through color object, since they has different ids 
        #assigned each time using color_dict_to_objects
        for key in data_colors:
            self.assertEqual(obs_data_colors[key].toHex(),\
                             data_colors[key].toHex())
        
        self.assertEqual(obs_data_color_order,self.data_color_order)
        
    def test_get_group_colors(self):
        """get_group_colors should iterate over color groups correctly."""

        data_colors = color_dict_to_objects(self.data_color_hsv)
        exp=(self.colors,data_colors,self.data_color_order)
        obs=get_group_colors(self.groups,self.colors,data_colors,\
                             self.data_color_order)
        
        self.assertEqual(obs,exp)

    def test_get_color(self):
        """get_color should get colors by several means"""
        sgc = lambda x: str(get_color(x))

        self.assertEqual(sgc('red1'), 'red1:#ff0000')
        self.assertEqual(sgc(('red1', '#00cc00')), 'red1:#00cc00')
        self.assertEqual(sgc(('red1', (120,100,50))), 'red1:#007f00')
        self.assertRaises(ValueError, sgc, 'xyz')
 
    def test_color_groups(self):
        """color_groups should iterate over color groups correctly."""
        data_colors = color_dict_to_objects(self.data_color_hsv)
       
        exp=None
        obs=color_groups(self.groups,data_colors,self.data_color_order)

        self.assertEqual(obs,exp)

    def test_make_color_dict(self):
        """make_color_dict: returns dict of named colors"""
        self.assertEqual(make_color_dict('red',(0,100,100),'white',(0,0,100),3),
            {   'redtowhite3_0':[0,100,100],
                'redtowhite3_1':[0,50,100],
                'redtowhite3_2':[0,0,100],
            })

    def test_process_colorby(self):
        """process_colorby: parses the cmd line and determines which columns \
from mapping file to color by"""
        self.colorby='Day'
        exp1={}
        exp1['Day']={'column':'Day'}
        obs1,obs2=process_colorby(self.colorby,self.data)

        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,self.data)

    def test_taxonomy_process_prefs(self):
        """taxonomy_process_prefs should return a taxonomy prefs file"""
        taxonomy_levels = [2,3,4]
        exp1={}
        exp1['level_2']={'column':'2'}
        exp1['level_3']={'column':3, 'colors':{'a':('red',(0, 100, 100))}}
        obs1=taxonomy_process_prefs(taxonomy_levels,exp1)
        
        exp2={}
        exp2['2']={'column':'2', 'colors':{}}
        exp2['3']={'column':'3', 'colors':{'a':('red',(0, 100, 100))}}
        exp2['4'] = {'column':'4','colors':{}}
        self.assertEqual(obs1,exp2)

        
#run tests if called from command line
if __name__ == "__main__":
    main()
