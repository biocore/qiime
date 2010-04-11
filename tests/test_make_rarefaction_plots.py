#!/usr/bin/env python
# File created on 10 Apr 2010
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from qiime.make_rarefaction_plots import make_plots,save_rarefaction_plots
from os.path import exists
from os import remove

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""

        self.data={}
        self.data['xaxis']=[10.0, 23.0, 36.0]
        self.data['yvals']={'Sample1': [1.3276140000000001, 2.1466460000000001], 'Sample2': [1.3276140000000001, 2.1466460000000001]}
        self.data['err']={'Sample1': [.1, .2], 'Sample2': [.2, .1]}
        self.xmax=140
        self.ymax=20
        self.ops=['Sample1','Sample2']
        self.mapping_category='SampleID'
        self.imagetype='png'
        self.resolution=70
        self.rtype='/tmp/SampleID.txt'
        self.data['map']=[['SampleID','Day'],['Sample1','Day1'],['Sample2',\
                          'Day1']]
        self.color_prefs={'SampleID': {'column': 'SampleID', 'color': {'Sample1': '#0000ff','Sample2': '#0000ff'}}}
        
        self.background_color='black'
        self.label_color='white'
        
        self.rares={'SampleID.txt': {'color': {'Sample1': '#0000ff','Sample2': '#0000ff'}, 'series': {'Sample1': [1.3276140000000001, 2.1466460000000001], 'Sample2': [1.3276140000000001, 2.1466460000000001]}, 'headers': ['', 'SampleID'], 'xaxis': [0,10.0], 'error': {'Sample1': [.1, .2], 'Sample2': [.2, .1]}}}
        self.fpath='/tmp/'
        self.output_dir='/tmp'
        
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
        
        self._paths_to_clean_up = []
    
    def tearDown(self):
        
        map(remove,self._paths_to_clean_up)
        
    def test_make_plots(self):
        """make_plots: creates HTML source for images"""

        filename1='/tmp/SampleID.png'
        self._paths_to_clean_up = [filename1]

        obs1 = make_plots(self.color_prefs, self.data, self.background_color, self.label_color, self.rares, self.ymax, self.output_dir, self.resolution, self.imagetype)

        self.assertEqual(obs1,exp_html)
        self.assertTrue(exists(filename1),'The png file was not created in \
#the appropriate location')


    
exp_html=\
'''\n<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n<html>\n<head>\n  <meta http-equiv="content-type" content="text/html;">\n  <title>Rarefaction Curves</title>\n<style type="text/css">\ntd.data{font-size:10px}\n</style>\n<script language="javascript" type="text/javascript">\n\nfunction changeCategory(SelObject){\n\nvar header_display=document.getElementById(\'rare_header\');\nheader_display.style.display=\'\';\n\nvar old_selected=document.getElementById(\'old_selected\');\n\ndata_display=document.getElementsByName(SelObject.value)\nfor (var i=0;i<data_display.length;i++){\n    data_display[i].style.display="";\n}\n\ndata_hide=document.getElementsByName(old_selected.value)\nfor (var i=0;i<data_hide.length;i++){\n    data_hide[i].style.display="none";\n}\n\nold_selected.value=SelObject.value;\n}\n</script>\n</head>\n<body>\n<form>\n<input id="old_selected" type="hidden" value="" \\>\n</form>\n<table><tr>\n<td><b>Select a Category:</b></td>\n<td>\n<select onchange="javascript:changeCategory(this)">\n<option></option>\n<option value="SampleID">SampleID</option>\n</select>\n</td>\n</table>\n<br>\n<table id="rare_plots">\n<tr name="SampleID" style="display: none;"><td colspan="20"><img src=".//SampleID.png" \\></td></tr>\n</table>\n<br>\n<table id="legend">\n<tr name="SampleID" style="display: none;"><td><b>SampleID Coloring:</b></td>\n<td bgcolor="#0000ff">Sample1</td>\n<td bgcolor="#00ff00">Sample2</td></tr>\n</table>\n<br>\n<table id="rare_data">\n<tr id="rare_header" style="display: none;">\n<td><b>Category</b></td>\n<td><b>SampleID</b></td>\n<td><b>Sequence Depth</b></td>\n<td><b>Rarefaction Ave.</b></td>\n<td><b>Error</b></td>\n</tr>\n<tr name="SampleID" style="display: none;"><td class="data" bgcolor="#0000ff">Sample1</td><td class="data">Sample1</td><td class="data">0</td><td class="data">1.327614</td><td class="data">0.1</td></tr>\n<tr name="SampleID" style="display: none;"><td class="data" bgcolor="#0000ff">Sample1</td><td class="data">Sample1</td><td class="data">1</td><td class="data">2.146646</td><td class="data">0.2</td></tr>\n<tr name="SampleID" style="display: none;"><td class="data" bgcolor="#00ff00">Sample2</td><td class="data">Sample2</td><td class="data">0</td><td class="data">1.327614</td><td class="data">0.2</td></tr>\n<tr name="SampleID" style="display: none;"><td class="data" bgcolor="#00ff00">Sample2</td><td class="data">Sample2</td><td class="data">1</td><td class="data">2.146646</td><td class="data">0.1</td></tr>\n</table>\n</body>\n</html>\n'''

if __name__ == "__main__":
    main()
