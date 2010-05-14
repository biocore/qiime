#!/usr/bin/env python
#file test_make_2d_plots.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"

import matplotlib
from matplotlib import use
use('Agg',warn=False)
from numpy import array
from os.path import exists
from StringIO import StringIO
from cogent.util.unit_test import TestCase, main
from os import remove
from qiime.make_2d_plots import (make_interactive_scatter,transform_xy_coords,
                                  draw_scatterplot,draw_pca_graph,
                                  extract_and_color_xy_coords,write_html_file,
                                  create_html_filename,
                                  convert_coord_data_to_dict,generate_xmap)
from qiime.colors import data_colors

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        
        self.props={"title": "PCA - P1 vs P2","ylabel":"P2","xlabel":"P1"}
        self.data={}
        self.data['coord']=[['Sample1','Sample2'],array([[-0.2,0.07],\
                            [-0.04,0.2]]),array([0.7,0.6]),\
                            array([25.00,30.00])]
        self.data['map']=[['#SampleID','Day'],['Sample1','Day1'],['Sample2',\
                          'Day1']]
        self.coord_tups = [("1", "2"), ("3", "2"), ("1", "3")]
        self.generate_eps=True
        self.data['alpha']=0.33
        self.groups={}
        self.groups['Day1']=['Sample1','Sample2']
        self.colors={}
        self.colors['Day1']='blue'
        self.prefs={}
        self.prefs['Sample']={}
        self.prefs['Sample']['column']='Day'
        self.data_color_hsv = {
              'aqua':     (180, 100, 100),
              'blue':     (240,100,100),
              'fuchsia':  (300,100,100),
              'gray':     (300,0,50.2),
              'green':    (120,100,50.2),
              'lime':     (120,100,100),
              'maroon':   (0,100,50.2),
              'olive':    (60,100,50.2),
              'purple':   (300,100,50.2),
              'red':      (0,100,100),
              'silver':   (0, 0, 75.3),
              'teal':     (180,100,50.2),
              'yellow':   (60,100,100),
        }
        self.data_color_order = ['blue','lime','red','aqua','fuchsia','yellow',\
                        'green','maroon','teal','purple','olive','silver','gray',[]]
        self.background_color='black'
        self.label_color='white'
        self.dir_path='/tmp/'
        self.data_file_link='/tmp/'
        self.xy_coords={}
        self.xy_coords['Sample1']=([-0.2], [0.07], ['Sample1: Day1'],\
                                   ['#0000ff'],['s'],[None],[None],[None])
        self.xy_coords['Sample2']=([-0.04], [0.2], ['Sample2: Day1'],\
                                   ['#0000ff'],['s'],[None],[None],[None])
                        
        self.coord_1='1'
        self.coord_2='2'
        
        self.p2d={}
        self.p2d['Sample1']=-0.2
        self.p2d['Sample2']=-0.04
        self.p1d={}
        self.p1d['Sample1']=0.07
        self.p1d['Sample2']=0.2
        self.all_cids={}
        self.all_cids=['Sample1: Day1', 'Sample2: Day1']
        self.all_xcoords=[100.79999999999998, 279.36000000000001] 
        self.all_ycoords=[54.000000000000014, 288.0]
        self.plot_label='SampleID'
        self.coords={'pc vector number':['Sample1','Sample2'],'1':\
                     array([-0.2,-0.04]),'2':array([0.07, 0.2])}
        self.x_len=4.5
        self.y_len=4.5
        self.size=20
        self.alpha=0.33
        self._paths_to_clean_up = []
    
    def tearDown(self):
        map(remove,self._paths_to_clean_up)

    def test_make_interactive_scatter(self):
        """make_interactive_scatter: creates HTML source for interactive \
images"""

        filename1='/tmp/P1_vs_P2_plot.png'
        filename2='/tmp/P1vsP2plot.eps.gz'

        self._paths_to_clean_up = [filename1,filename2]

        obs1,obs2,obs3=make_interactive_scatter(self.plot_label,self.dir_path,
                                self.data_file_link,self.background_color,
                                self.label_color,None,self.alpha, 
                                self.xy_coords,self.props, 
                                self.x_len, self.y_len, self.size,
                                draw_axes=False, generate_eps=True)

        self.assertEqual(obs1,expsrcmap1)
        self.assertEqual(obs2,expimgmap1)
        self.assertEqual(obs3,expeps1)
        self.assertTrue(exists(filename1),'The png file was not created in \
the appropriate location')
        self.assertTrue(exists(filename2),'The eps file was not created in \
the appropriate location')
    
    def test_generate_xmap(self):
        """generate_xmap: generates the html area map"""
        exp2=360
        exp3=360
        obs1,obs2,obs3=generate_xmap(self.x_len,self.y_len,self.all_cids,\
                                     self.all_xcoords,self.all_ycoords)
        self.assertEqual(obs1,exparea)
        self.assertEqual(obs2,exp2)
        self.assertEqual(obs3,exp3)

    def test_draw_scatterplot(self):
        """draw_scatterplot: draws the matplotlib scatterplot"""
        exp=array([[-0.04, 0.2 ]])

        sc_plot = draw_scatterplot(self.props,self.xy_coords,self.x_len,\
                                   self.y_len,self.size,
                                   self.background_color,self.label_color, None,
                                   self.alpha)
        obs=sc_plot.get_offsets()

        self.assertEqual(obs,exp)
        
    def test_transform_xy_coords(self):
        """transform_xy_coords: transforms the xy coords from the matplotlib \
plot into html spatial coords which allows for mouseovers"""
        sc_plot = draw_scatterplot(self.props,self.xy_coords,self.x_len,\
                                   self.y_len, self.size,
                                   self.background_color,self.label_color, None,
                                   self.alpha)
                               
        obs1,obs2,obs3=transform_xy_coords(self.xy_coords,sc_plot)
        
        self.assertEqual(obs1,self.all_cids)
        self.assertEqual(obs2,self.all_xcoords)
        self.assertEqual(obs3,self.all_ycoords)
        
    def test_draw_pca_graph(self):
        """draw_pca_graph: draws the matplotlib figure"""

        filename1='/tmp/P1_vs_P2_plot.png'
        filename2='/tmp/P1vsP2plot.eps.gz'

        self._paths_to_clean_up = [filename1,filename2]
        
        obs1,obs2=draw_pca_graph(self.plot_label,self.dir_path,
                                 self.data_file_link,self.coord_1,self.coord_2,
                                 None,None,None,None,
                                 self.data,self.prefs,self.groups,self.colors,
                                 self.background_color,self.label_color,
                                 data_colors,self.data_color_order,
                                 generate_eps=True)
                                 

        self.assertEqual(obs1,expsrcmap2+expimgmap2)
        self.assertEqual(obs2,expeps2)
        self.assertTrue(exists(filename1),'The png file was not created in \
the appropriate location')
        self.assertTrue(exists(filename2),'The eps file was not created in \
the appropriate location')

    def test_extract_and_color_xy_coords(self):
        """extract_and_color_xy_coords: gets coords from coords file and \
associates colors to those coords based on its group"""
        
        obs=extract_and_color_xy_coords(self.p1d,self.p2d,None,None,None,self.colors,
                                        data_colors,self.groups,self.coords)
        
        self.assertFloatEqual(obs,self.xy_coords)
        
    def test_create_html_filename(self):
        """create_html_filename: using the pca filename, generates an html \
filename for the plots"""
        
        exp='test_2D.html'
        obs=create_html_filename(coord_filename='test',name_ending='_2D.html')
        self.assertEqual(obs,exp)
        
    def test_convert_coord_data_to_dict(self):
        """convert_coord_data_to_dict: converts the coords list into a \
dictionary"""
        
        exp1={'pc vector number':['Sample1','Sample2'],'1':array([-0.2,-0.04]),\
                '2':array([0.07, 0.2])}
        exp2={'1':[25.00],'2':[30.00],}
        obs1,obs2=convert_coord_data_to_dict(self.data)
        
        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,exp2)
        
    def test_write_html_file(self):
        "Write html and make sure it gets cleaned up"""
        filename1='/tmp/test.html'
        
        self._paths_to_clean_up = [filename1]
        
        write_html_file('Test','/tmp/test.html')
        
        self.assertTrue(exists(filename1),'The file was not created in \
the appropriate location')

#expected results for the unit testing       
exparea=['<AREA shape="circle" coords="100,306,5" href="#Sample1: Day1"  onmouseover="return overlib(\'Sample1: Day1\');" onmouseout="return nd();">\n', '<AREA shape="circle" coords="279,72,5" href="#Sample2: Day1"  onmouseover="return overlib(\'Sample2: Day1\');" onmouseout="return nd();">\n']

expsrcmap1 = '<img src="/tmp/P1_vs_P2_plot.png" border="0" ismap usemap="#pointsSampleID12" width="360" height="360" />\n'
expimgmap1 = '\n<MAP name="pointsSampleID12">\n\
<AREA shape="circle" coords="100,306,5" href="#Sample1: Day1"  onmouseover="return overlib(\'Sample1: Day1\');" onmouseout="return nd();">\n\
<AREA shape="circle" coords="279,72,5" href="#Sample2: Day1"  onmouseover="return overlib(\'Sample2: Day1\');" onmouseout="return nd();">\n\n\
</MAP>\n'
expeps1='<a href="/tmp/P1vsP2plot.eps.gz" >Download Figure</a>'

expsrcmap2 = '<img src="/tmp/P1_vs_P2_plot.png" border="0" ismap usemap="#pointsSampleID12" width="360" height="360" />\n'
expimgmap2 = '\n<MAP name="pointsSampleID12">\n\
<AREA shape="circle" coords="100,243,5" href="#Sample1: Day1"  onmouseover="return overlib(\'Sample1: Day1\');" onmouseout="return nd();">\n\
<AREA shape="circle" coords="279,93,5" href="#Sample2: Day1"  onmouseover="return overlib(\'Sample2: Day1\');" onmouseout="return nd();">\n\n\
</MAP>\n'
expeps2='<a href="/tmp/P1vsP2plot.eps.gz" >Download Figure</a>'

#run tests if called from command line
if __name__ == "__main__":
    main()
