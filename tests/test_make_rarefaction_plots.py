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
from qiime.make_rarefaction_plots import make_plots, \
                    save_ave_rarefaction_plots,save_single_rarefaction_plots, \
                    get_rarefaction_data,make_html,make_averages, \
                    save_rarefaction_data,ave_seqs_per_sample, \
                    make_error_series
from os.path import exists
from os import remove
from qiime.parse import parse_rarefaction, parse_mapping_file
from qiime.colors import color_dict_to_objects
from shutil import rmtree
from qiime.util import create_dir

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""

        self.data={}
        self.data['xaxis']=[10.0, 23.0]
        self.sample_dict={'Sample1':{0.00: [1.3276140000000001], \
                                     10.0:[2.1466460000000001]}}
        self.data['yvals']={'Sample1': [1.3276140000000001, 2.1466460000000001],
                            'Sample2': [1.3276140000000001, 2.1466460000000001]}
        self.data['err']={'Sample1': [.1, .2], 'Sample2': [.2, .1]}
        self.xmax=140
        self.ymax=20
        self.ops=['Sample1','Sample2']
        self.mapping_category='SampleID'
        self.imagetype='png'
        self.resolution=70
        self.data['map']=[['SampleID','Day'],['Sample1','Day1'],['Sample2',\
                          'Day1']]
        self.color_prefs={'SampleID': {'column': 'SampleID', 'color': \
                          {'Sample1': '#0000ff','Sample2': '#0000ff'}}}
        self.groups={'Sample1':['Sample1']}
        self.background_color='black'
        self.label_color='white'
        self.labelname='SampleID'
        self.rare_data={'color': {'Sample1': '#00ff7f'}, \
            'series': {'Sample1': [2.0515300000000001, 2.8053979999999998],}, \
             'headers': ['test.txt','SampleID'], 'xaxis': [10.0, 23.0], \
             'error': {'Sample1': [0.0, 0.0]}, 'options': ['Sample1']}
        self.fpath='/tmp/'
        self.output_dir='/tmp/'
        self.metric_name='test'
        self._paths_to_clean_up = []
        self._folders_to_cleanup = []
        self.rarefaction_file_data=[[10.0, 0.0, 1.0], [10.0, 1.0, 3.0]]
        d = {'redtowhite3_0':[0,100,100],'redtowhite3_1':[0,50,100]}
        self.data_colors = color_dict_to_objects(d)
        self.colors={'Sample1':'redtowhite3_0','Sample2':'redtowhite3_1'}
        self.colors2={'Sample1':'redtowhite3_0'}
        self.mappingfile = ['#SampleID\tSex\tAge',
                            '123\tF\t32',
                            '234\tM\t30',
                            '345\tM\t32']
        #self.p_mappingfile = parse_mapping_file(self.mappingfile,\
        #                                            strip_quotes=True)
        self.rarefactionfile=[\
                    '\tsequences per sample\titeration\t123\t234\t345',
                    'rare10.txt\t10\t0\t1.99181\t0.42877\t2.13996',
                    'rare10.txt\t10\t1\t2.07163\t0.42877\t2.37055',
                    'rare310.txt\t310\t0\t8.83115\t0.42877\t11.00725',
                    'rare310.txt\t310\t1\t10.05242\t0.42877\t8.24474',
                    'rare610.txt\t610\t0\t12.03067\t0.42877\t11.58928',
                    'rare610.txt\t610\t1\t12.9862\t0.42877\t11.58642']
                    
        self.rares = {'test.txt': (['', 'sequences per sample', 'iteration', \
                      'Sample1', 'Sample2'], [], ['rare1.txt', 'rare2.txt'], \
                      [[10.0, 2.0, 7.0, 7.0, 9.0], [10.0, 2.0, 7.0, 7.0, 9.0]])}
        self.col_headers, self.comments, self.rarefaction_fns, \
        self.rarefaction_data = parse_rarefaction(self.rarefactionfile)
        self.matrix, self.seqs_per_samp, self.sampleIDs = \
        get_rarefaction_data(self.rarefaction_data, self.col_headers)
        self.ave_seqs_per_sample = {'123':[2.03172,9.4417849999999994,\
        12.508435],'234':[0.42876999999999998,0.42876999999999998,\
        0.42876999999999998],'345':[2.255255,9.625995,11.58785]}
        self.collapsed_ser_sex = {'M':[1.3420125000000001,5.0273824999999999,\
        6.0083099999999998], 'F':[2.03172,9.4417849999999994,12.508435]}
        self.err_ser_sex = {'M':[0.91324250000000007,4.5986124999999998,\
        5.5795399999999997],'F':[0.0,0.0,0.0]}
        self.rarefaction_legend_mat_init={'test': {'SampleID': {}}}
        self.col_headers2=['', 'sequences per sample', 'iteration', 'Sample1', \
                           'Sample2']
        self.rarefaction_legend_mat_init_single_rare={'test': {'SampleID': \
            {'Sample1': {'link': 'all_other_plots/testSampleIDSample1.png',\
            'groupcolor': '#0000ff', 'groupsamples': {'Sample1': '#0000ff'}}}}}
        self.rarefaction_legend_mat_returned={'test': {'SampleID': {'Sample1':\
            {'link': 'tmp/testSampleID.png',  'groupcolor': '#00ff7f', \
            'groupsamples': {'Sample1': '#007fff'}},'link': \
            'tmp/testSampleID.png'}}}
        self.rarefaction_data_mat={'SampleID': {'Sample1': {'test': \
                                   {'ave': ['     1.361', '     2.175'],'err':\
                                   ['     0.000', '     0.000']}}}}
    
    def tearDown(self):
        '''This function removes the generated files'''
        map(remove,self._paths_to_clean_up)
        map(rmtree,self._folders_to_cleanup)
    
    def test_save_ave_rarefaction_plots(self):
        '''save_ave_rarefaction_plots: this tests the functionality of 
           creating a rarefaction plot based on the average rarefaction 
           values'''
           
        filename1='/tmp/test1SampleID.png'
        self._paths_to_clean_up = [filename1]

        exp={'test': {'SampleID': {'link': 'average_plots/testSampleID.png'}}}
        obs=save_ave_rarefaction_plots(self.data['xaxis'], \
                    self.data['yvals'],self.data['err'], \
                    self.xmax, self.ymax, self.ops, self.mapping_category, \
                    self.imagetype, self.resolution, self.data_colors, \
                    self.colors, '/tmp/test1',self.background_color, \
                    self.label_color,self.rarefaction_legend_mat_init, \
                    self.metric_name)
                    
        self.assertEqual(obs,exp)
        self.assertTrue(exists(filename1))
                    
    def test_get_rarefaction_data(self):
        '''get_rarefaction_data: This tests the functionality of taking a 
           rarefaction file and returns three values'''
            
        obs1,obs2,obs3 = get_rarefaction_data(self.rarefaction_file_data, \
                                                self.col_headers2)
                                                    
        exp1=[[1.0, 3.0]]
        exp2=[10.0, 10.0]
        exp3=['Sample1', 'Sample2']
        
        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,exp2)
        self.assertEqual(obs3,exp3)
        
    def test_make_html(self):
        '''make_html: this tests the output html format'''
        obs=make_html(self.rarefaction_legend_mat_returned, \
                         self.rarefaction_data_mat,self.data['xaxis'], \
                         self.imagetype)
                         
        self.assertEqual(obs,exp_html)
    
    def test_make_averages(self):
        '''make_averages: this tests the main function in make_rarefaction_plots
           and returns an html output'''
           
        filename1='/tmp/all_other_plots/testSampleIDSample1.png'
        filename2='/tmp/all_other_plots/testSampleIDSample2.png'
        folder1='/tmp/all_other_plots'
        filename3='/tmp/average_tables/testSampleID.txt'
        folder2='/tmp/average_tables'
        filename4='/tmp/average_plots/testSampleID.png'
        folder3='/tmp/average_plots'
        self._paths_to_clean_up = [filename1,filename2,filename3,filename4]
        self._folders_to_cleanup = [folder1,folder2,folder3]
        
        obs=make_averages(self.color_prefs,self.data,self.background_color, \
                          self.label_color,self.rares,self.output_dir, \
                          self.resolution,self.imagetype,None)
                          
        self.assertEqual(obs,exp_make_avg)
        self.assertTrue(exists(filename1))
        self.assertTrue(exists(filename2))
        self.assertTrue(exists(filename3))
        self.assertTrue(exists(filename4))
        self.assertTrue(exists(folder1)) 
        self.assertTrue(exists(folder2)) 
        self.assertTrue(exists(folder3)) 
        
    def test_save_rarefaction_data(self):
        '''save_rarefaction_data: This tests the rarefaction average output'''
        
        exp=['# test.txt\n', '# SampleID\n', 'xaxis: 10.0\t23.0\t\n', \
        'xmax: 140\n', '>> Sample1\n', 'color #ff0000\n', 'series ', \
        'nan\n', 'error ', '0.0\n']
                     
        obs=save_rarefaction_data(self.ave_seqs_per_sample,self.data['xaxis'],\
                                  self.xmax,self.mapping_category,self.colors2,\
                                  'test.txt',self.data_colors,self.groups)
                    
        self.assertEqual(obs,exp)
        
    def test_save_single_rarefaction_plots(self):
        '''save_single_rarefaction_plots: this generates a plot with raw
           rarefaction data on a plot and tests whether a file is generated'''
           
        filename1='/tmp/testSampleID.png'
        self._paths_to_clean_up = [filename1]

        exp={'test': {'SampleID': {'Sample1': {'groupcolor': '#0000ff', 'link':\
                     'all_other_plots/testSampleIDSample1.png', 'groupsamples':\
                     {'Sample1': '#ff0000'}}}}}
                     
        obs=save_single_rarefaction_plots(['Sample1'],self.sample_dict, \
                    self.imagetype,self.metric_name,'SampleID', \
                    self.data_colors, self.colors, '/tmp/testSampleID', \
                    self.background_color, self.label_color, self.resolution, \
                    self.ymax,self.xmax,'Sample1', \
                    self.rarefaction_legend_mat_init_single_rare)
                    
        self.assertEqual(obs,exp)
        self.assertTrue(exists(filename1))
    
    def test_ave_seqs_per_sample(self):
        '''ave_seqs_per_sample: this tests getting the average seqs per 
           sample'''
           
        test = ave_seqs_per_sample(self.matrix,self.seqs_per_samp,\
        self.sampleIDs)
        
        self.assertEqual(test, self.ave_seqs_per_sample)
    
    def test_make_error_series(self):
        '''make_error_series: this tests whether the errors were correctly
           calculated'''
           
        groups={'M': ['234', '345'], 'F': ['123']}
        
        test = make_error_series(self.ave_seqs_per_sample, groups)
        
        self.assertEqual(test[0], self.collapsed_ser_sex)
        self.assertEqual(test[1], self.err_ser_sex)

    def test_make_plots(self):
        """make_plots: tests whether the average plots are generated and if
           dictionary for the html generation is properly formatted"""

        filename1='/tmp/test/testSampleID.png'
        folder1='/tmp/test/'
        
        self._paths_to_clean_up = [filename1]
        self._folders_to_cleanup=[folder1]
        
        exp1={'SampleID': {'Sample1': {'test': {'ave': ['     1.361', \
        '     2.175', '     2.052', '     2.805'], 'err': ['     0.000', \
        '     0.000', '     0.000', '     0.000']}}}}
        
        exp2={'test': {'SampleID': {'link': 'average_plots/testSampleID.png', \
        'Sample1': {'groupcolor': '#00ff7f', 'link': 'tmp/testSampleID.png', 
        'groupsamples': {'Sample1': '#007fff'}}}}}
        
        create_dir('/tmp/test/',False)
        
        obs1,obs2 = make_plots(self.background_color,self.label_color, \
                          self.rare_data,self.ymax, self.xmax,'/tmp/test/', \
                          self.resolution, self.imagetype,self.groups,\
                          self.colors,self.data_colors,self.metric_name,\
                          self.labelname,self.rarefaction_data_mat, \
                          self.rarefaction_legend_mat_returned)
        
        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,exp2)
        self.assertTrue(exists(filename1))
        self.assertTrue(exists(folder1))
        
    '''   
    #These functions are currently not being used
    
    def test_get_overall_averages(self):
        test = get_overall_averages(self.ave_seqs_per_sample,self.sampleIDs)
        self.assertEqual(test, self.overall_averages)
        
    def test_is_max_category_ops_neg(self):
            test = is_max_category_ops(self.p_mappingfile, 'Sex')
            self.assertEqual(test[0], False)

    def test_is_max_category_ops_pos(self):
            test = is_max_category_ops(self.p_mappingfile, 'SampleID')
            self.assertEqual(test[0], True)
    '''

    
exp_html=\
'''\n<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n<html>\n<head>\n  <meta http-equiv="content-type" content="text/html;">\n  <title>Rarefaction Curves</title>\n<style type="text/css">\ntd.data{font-size:10px}\ntd.headers{font-size:12px}\ntable{border-spacing:0px}\n.removed{display:none;}\n.expands{cursor:pointer; cursor:hand;}\n.child1 td:first-child{padding-left: 3px;}\n</style>\n<script language="javascript" type="text/javascript">\n\nfunction change_plot(object){\nplot_img=document.getElementById(\'plots\');\nplot_img.src=object\n}\n\nfunction changeMetric(SelObject){\n    var old_category=document.getElementById(\'category\');\n    var old_metric=document.getElementById(\'metric\');\n    var imagetype=document.getElementById(\'imagetype\').value;\n    var legend=document.getElementById(\'legend\');\n\n    if (old_category.value != \'\'){\n        plot_img=document.getElementById(\'plots\');\n        plot_loc=\'./average_plots/\'+SelObject.value+old_category.value+imagetype\n        plot_img.src=plot_loc\n        plot_img.style.display="";\n        legend.style.display="";\n         data_display=document.getElementsByName(SelObject.value+old_category.value)\n        for (var i=0;i<data_display.length;i++){\n            data_display[i].style.display="";\n        }\n        data_hide=document.getElementsByName(old_metric.value+old_category.value)\n        for (var i=0;i<data_hide.length;i++){\n            data_hide[i].style.display="none";\n        }\n        \n        data_display=document.getElementsByName(old_category.value)\n        for (var i=0;i<data_display.length;i++){\n            data_display[i].style.display="";\n        }\n    }\n    \nold_metric.value=SelObject.value;\n}\n\nfunction changeCategory(SelObject){\n    var old_category=document.getElementById(\'category\');\n    var old_metric=document.getElementById(\'metric\');\n    var imagetype=document.getElementById(\'imagetype\').value;\n    var legend=document.getElementById(\'legend\');\n\n    if (old_metric.value != \'\'){\n        plot_img=document.getElementById(\'plots\');\n        plot_loc=\'./average_plots/\'+old_metric.value+SelObject.value+imagetype;\n        plot_img.src=plot_loc;\n        plot_img.style.display="";\n        legend.style.display="";\n        \n        data_display=document.getElementsByName(SelObject.value)\n        for (var i=0;i<data_display.length;i++){\n            data_display[i].style.display="";\n        }\n\n        data_hide=document.getElementsByName(old_category.value)\n        for (var i=0;i<data_hide.length;i++){\n            data_hide[i].style.display="none";\n        }\n\n        data_display=document.getElementsByName(old_metric.value+SelObject.value)\n        for (var i=0;i<data_display.length;i++){\n            data_display[i].style.display="";\n        }\n\n        data_hide=document.getElementsByName(old_metric.value+old_category.value)\n        for (var i=0;i<data_hide.length;i++){\n            data_hide[i].style.display="none";\n        }\n    }\nold_category.value=SelObject.value;\n}\n\nfunction toggle(){\n    if (arguments[0].innerHTML==\'\\u25B6\'){\n        arguments[0].innerHTML=\'\\u25BC\'\n    }else{\n        arguments[0].innerHTML=\'\\u25B6\'\n    }\n    for(var i=1; i<arguments.length; i++){\n        with(document.getElementById(arguments[i])){\n            if(className.indexOf(\'removed\') > -1){\n                className = className.replace(\'removed\');\n            }else{\n                className += \' removed\';\n            }\n        }\n    }\n}\n</script>\n\n\n\n</head>\n<body>\n<form>\n<input id="metric" type="hidden"></input>\n<input id="category" type="hidden"></input>\n<input id="imagetype" type="hidden" value=".png"></input>\n</form>\n<table><tr>\n<td><b>Select a Metric:</b></td>\n<td>\n<select onchange="javascript:changeMetric(this)">\n<option></option>\n<option value="test">test</option>\n</select>\n</td>\n<td><b>Select a Category:</b></td>\n<td>\n<select onchange="javascript:changeCategory(this)">\n<option></option>\n<option value="SampleID">SampleID</option>\n</select>\n</td>\n</table>\n<br>\n<table><tr><td><img id="plots" style="display: none; height: 400px;" \\></td><td id="legend" style="display: none;"><b>Legend</b><div STYLE="border: thin black solid; height: 300px; width: 150px; font-size: 12px; overflow: auto;"><table><tr id="testSampleID" class="expands" onclick="toggle(testSampleIDSample1,\'testSampleIDSample1\')" name="testSampleID" style="display: none;" onmouseover="javascript:change_plot(\'./tmp/testSampleID.png\')" onmouseout="javascript:change_plot(\'./tmp/testSampleID.png\')"><td id="testSampleIDSample1"  class="data">&#x25B6;</td><td class="data" bgcolor="#00ff7f"><b>Sample1</b></td></tr>\n<tr id="testSampleIDSample1" name="testSampleID" class="removed child1"><td class="data" align="right">&#x221F;</td><td bgcolor="#007fff" class="data" style="border:thin white solid;"><b>Sample1</b></td></tr></table></div></td></tr></table>\n<br>\n<table id="rare_data">\n<tr name="SampleID" style="display: none;"><td class="headers">SampleID</td><td class="headers">Seqs/Sample</td>\n<td class="headers">test Ave.</td><td class="headers">test Err.</td>\n</tr>\n<tr name="SampleID" style="display: none;"></tr>\n<tr name="SampleID" style="display: none;">\n<td class="data" bgcolor="#00ff7f">Sample1</td><td class="data">10.0</td>\n<td class="data">     1.361</td><td class="data">     0.000</td>\n<tr name="SampleID" style="display: none;">\n<td class="data" bgcolor="#00ff7f">Sample1</td><td class="data">23.0</td>\n<td class="data">     2.175</td><td class="data">     0.000</td>\n</tr>\n</table>\n</body>\n</html>\n'''

exp_make_avg=\
'''\n<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n<html>\n<head>\n  <meta http-equiv="content-type" content="text/html;">\n  <title>Rarefaction Curves</title>\n<style type="text/css">\ntd.data{font-size:10px}\ntd.headers{font-size:12px}\ntable{border-spacing:0px}\n.removed{display:none;}\n.expands{cursor:pointer; cursor:hand;}\n.child1 td:first-child{padding-left: 3px;}\n</style>\n<script language="javascript" type="text/javascript">\n\nfunction change_plot(object){\nplot_img=document.getElementById(\'plots\');\nplot_img.src=object\n}\n\nfunction changeMetric(SelObject){\n    var old_category=document.getElementById(\'category\');\n    var old_metric=document.getElementById(\'metric\');\n    var imagetype=document.getElementById(\'imagetype\').value;\n    var legend=document.getElementById(\'legend\');\n\n    if (old_category.value != \'\'){\n        plot_img=document.getElementById(\'plots\');\n        plot_loc=\'./average_plots/\'+SelObject.value+old_category.value+imagetype\n        plot_img.src=plot_loc\n        plot_img.style.display="";\n        legend.style.display="";\n         data_display=document.getElementsByName(SelObject.value+old_category.value)\n        for (var i=0;i<data_display.length;i++){\n            data_display[i].style.display="";\n        }\n        data_hide=document.getElementsByName(old_metric.value+old_category.value)\n        for (var i=0;i<data_hide.length;i++){\n            data_hide[i].style.display="none";\n        }\n        \n        data_display=document.getElementsByName(old_category.value)\n        for (var i=0;i<data_display.length;i++){\n            data_display[i].style.display="";\n        }\n    }\n    \nold_metric.value=SelObject.value;\n}\n\nfunction changeCategory(SelObject){\n    var old_category=document.getElementById(\'category\');\n    var old_metric=document.getElementById(\'metric\');\n    var imagetype=document.getElementById(\'imagetype\').value;\n    var legend=document.getElementById(\'legend\');\n\n    if (old_metric.value != \'\'){\n        plot_img=document.getElementById(\'plots\');\n        plot_loc=\'./average_plots/\'+old_metric.value+SelObject.value+imagetype;\n        plot_img.src=plot_loc;\n        plot_img.style.display="";\n        legend.style.display="";\n        \n        data_display=document.getElementsByName(SelObject.value)\n        for (var i=0;i<data_display.length;i++){\n            data_display[i].style.display="";\n        }\n\n        data_hide=document.getElementsByName(old_category.value)\n        for (var i=0;i<data_hide.length;i++){\n            data_hide[i].style.display="none";\n        }\n\n        data_display=document.getElementsByName(old_metric.value+SelObject.value)\n        for (var i=0;i<data_display.length;i++){\n            data_display[i].style.display="";\n        }\n\n        data_hide=document.getElementsByName(old_metric.value+old_category.value)\n        for (var i=0;i<data_hide.length;i++){\n            data_hide[i].style.display="none";\n        }\n    }\nold_category.value=SelObject.value;\n}\n\nfunction toggle(){\n    if (arguments[0].innerHTML==\'\\u25B6\'){\n        arguments[0].innerHTML=\'\\u25BC\'\n    }else{\n        arguments[0].innerHTML=\'\\u25B6\'\n    }\n    for(var i=1; i<arguments.length; i++){\n        with(document.getElementById(arguments[i])){\n            if(className.indexOf(\'removed\') > -1){\n                className = className.replace(\'removed\');\n            }else{\n                className += \' removed\';\n            }\n        }\n    }\n}\n</script>\n\n\n\n</head>\n<body>\n<form>\n<input id="metric" type="hidden"></input>\n<input id="category" type="hidden"></input>\n<input id="imagetype" type="hidden" value=".png"></input>\n</form>\n<table><tr>\n<td><b>Select a Metric:</b></td>\n<td>\n<select onchange="javascript:changeMetric(this)">\n<option></option>\n<option value="test">test</option>\n</select>\n</td>\n<td><b>Select a Category:</b></td>\n<td>\n<select onchange="javascript:changeCategory(this)">\n<option></option>\n<option value="SampleID">SampleID</option>\n</select>\n</td>\n</table>\n<br>\n<table><tr><td><img id="plots" style="display: none; height: 400px;" \\></td><td id="legend" style="display: none;"><b>Legend</b><div STYLE="border: thin black solid; height: 300px; width: 150px; font-size: 12px; overflow: auto;"><table><tr id="testSampleID" class="expands" onclick="toggle(testSampleIDSample1,\'testSampleIDSample1\')" name="testSampleID" style="display: none;" onmouseover="javascript:change_plot(\'./all_other_plots/testSampleIDSample1.png\')" onmouseout="javascript:change_plot(\'./average_plots/testSampleID.png\')"><td id="testSampleIDSample1"  class="data">&#x25B6;</td><td class="data" bgcolor="#0000ff"><b>Sample1</b></td></tr>\n<tr id="testSampleIDSample1" name="testSampleID" class="removed child1"><td class="data" align="right">&#x221F;</td><td bgcolor="#0000ff" class="data" style="border:thin white solid;"><b>Sample1</b></td></tr>\n<tr id="testSampleID" class="expands" onclick="toggle(testSampleIDSample2,\'testSampleIDSample2\')" name="testSampleID" style="display: none;" onmouseover="javascript:change_plot(\'./all_other_plots/testSampleIDSample2.png\')" onmouseout="javascript:change_plot(\'./average_plots/testSampleID.png\')"><td id="testSampleIDSample2"  class="data">&#x25B6;</td><td class="data" bgcolor="#00ff00"><b>Sample2</b></td></tr>\n<tr id="testSampleIDSample2" name="testSampleID" class="removed child1"><td class="data" align="right">&#x221F;</td><td bgcolor="#00ff00" class="data" style="border:thin white solid;"><b>Sample2</b></td></tr></table></div></td></tr></table>\n<br>\n<table id="rare_data">\n<tr name="SampleID" style="display: none;"><td class="headers">SampleID</td><td class="headers">Seqs/Sample</td>\n<td class="headers">test Ave.</td><td class="headers">test Err.</td>\n</tr>\n<tr name="SampleID" style="display: none;"></tr>\n<tr name="SampleID" style="display: none;">\n<td class="data" bgcolor="#0000ff">Sample1</td><td class="data">10.0</td>\n<td class="data">     7.000</td><td class="data">     0.000</td>\n<tr name="SampleID" style="display: none;">\n<td class="data" bgcolor="#00ff00">Sample2</td><td class="data">10.0</td>\n<td class="data">     7.000</td><td class="data">     0.000</td>\n</tr>\n</table>\n</body>\n</html>\n'''

if __name__ == "__main__":
    main()
