#!/usr/bin/env python
# File created on 10 Apr 2010
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"
 

from cogent.util.unit_test import TestCase, main
from qiime.make_rarefaction_plots import make_plots, \
                    save_single_ave_rarefaction_plots,save_single_rarefaction_plots, \
                    get_rarefaction_data,make_html,make_averages, \
                    save_rarefaction_data,ave_seqs_per_sample, \
                    make_error_series,save_ave_rarefaction_plots
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
        self.data['xaxis']=[10.0]
        self.sample_dict={'Sample1':{10.00: [1.3276140000000001]}}
        self.data['yvals']={'Sample1': [1.3276140000000001]}
        self.data['err']={'Sample1': [.1]}
        self.xmax=140
        self.ymax=20
        self.std_type='stddev'
        self.ops=['Sample1']
        self.mapping_category='SampleID'
        self.imagetype='png'
        self.resolution=70
        self.mapping_lookup={'SampleID-Sample1':'col_0_row_0'}
        self.data['map']=[['SampleID','Day'],['Sample1','Day1']]
        self.color_prefs={'SampleID': {'column': 'SampleID', 'color': \
                          {'Sample1': '#ff0000'}}}
        self.groups={'Sample1':['Sample1']}
        self.background_color='black'
        self.label_color='white'
        self.labelname='SampleID'
        self.rare_data={'color': {'Sample1': '#ff0000'}, \
            'series': {'Sample1': [2.0515300000000001],}, \
             'headers': ['test.txt','SampleID'], 'xaxis': [10.0], \
             'error': {'Sample1': [0.0]}, 'options': ['Sample1']}
        self.fpath='/tmp/'
        self.output_dir='/tmp/'
        self.metric_name='test'
        self._paths_to_clean_up = []
        self._folders_to_cleanup = []
        self.rarefaction_file_data=[[10.0, 0.0, 1.0], [10.0, 1.0, 3.0]]
        d = {'redtowhite3_0':'#7fff00','redtowhite3_1':'#7fff00'}
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
                      'Sample1'], [], ['rare1.txt', 'rare2.txt'], \
                      [[10.0, 2.0, 7.0, 7.0, 9.0], [10.0, 2.0, 7.0, 7.0, 9.0]])}
        self.col_headers, self.comments, self.rarefaction_fns, \
        self.rarefaction_data = parse_rarefaction(self.rarefactionfile)
        self.matrix, self.seqs_per_samp, self.sampleIDs = \
        get_rarefaction_data(self.rarefaction_data, self.col_headers)
        self.ave_seqs_per_sample1 = {'Sample1':[2.03172,9.4417849999999994,\
        12.508435]}
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
                           
        self.rarefaction_data_mat={'SampleID': {'Sample1': {'test': {'ave': ['     7.000'], 'err': ['       nan']}}}}
       
        self.rarefaction_legend_mat={'test': {'samples': {'Sample1': {'color': '#ff0000', 'link': 'html_plots/testcol_0_row_0.png'}}, 'groups': {'SampleID': {'Sample1': {'groupcolor': '#ff0000', 'groupsamples': ['Sample1']}}}}}
        self.exp_err_series_ave={'M': [1.571915, 6.49885, 8.1750183333333339]}
    
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

        exp=None
        
        obs=save_ave_rarefaction_plots(self.data['xaxis'], \
                    self.data['yvals'],self.data['err'], \
                    self.xmax, self.ymax, self.ops, self.mapping_category, \
                    self.imagetype, self.resolution, self.data_colors, \
                    self.colors, '/tmp/test1',self.background_color, \
                    self.label_color,self.metric_name)
                
        self.assertEqual(obs,exp)
        self.assertTrue(exists(filename1))
                    
    def test_save_single_ave_rarefaction_plots(self):
        '''save_single_ave_rarefaction_plots: this tests the functionality of 
           creating a rarefaction plot based on the average rarefaction 
           values'''
           
        filename1='/tmp/test1col_0_row_0_ave.png'
        self._paths_to_clean_up = [filename1]

        exp={'test': {'groups': {'SampleID': {'Sample1': {'groupcolor': '#ff0000', 'ave_link': 'html_plots/testcol_0_row_0_ave.png', 'groupsamples': ['Sample1']}}}, 'samples': {'Sample1': {'color': '#ff0000', 'link': 'html_plots/testcol_0_row_0.png'}}}}
        
        obs=save_single_ave_rarefaction_plots(self.data['xaxis'], \
                    self.data['yvals'],self.data['err'], \
                    self.xmax, self.ymax, self.ops, self.mapping_category, \
                    self.imagetype, self.resolution, self.data_colors, \
                    self.colors, '/tmp/test1',self.background_color, \
                    self.label_color,self.rarefaction_legend_mat, \
                    self.metric_name,self.mapping_lookup)
                    
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
        obs=make_html(self.rarefaction_legend_mat, \
                         self.rarefaction_data_mat,self.data['xaxis'], \
                         self.imagetype,self.mapping_lookup)
                         
        self.assertEqual(obs,exp_html)
    
    def test_make_averages(self):
        '''make_averages: this tests the main function in make_rarefaction_plots
           and returns an html output'''
           
        filename1='/tmp/html_plots/testcol_0_row_0_ave.png'
        folder1='/tmp/html_plots'
        filename3='/tmp/average_tables/testSampleID.txt'
        folder2='/tmp/average_tables'
        filename4='/tmp/average_plots/testSampleID.png'
        folder3='/tmp/average_plots'
        self._paths_to_clean_up = [filename1,filename3,filename4]
        self._folders_to_cleanup = [folder1,folder2,folder3]
        
        obs=make_averages(self.color_prefs,self.data,self.background_color, \
                          self.label_color,self.rares,self.output_dir, \
                          self.resolution,self.imagetype,None,False, \
                          self.std_type)
                          
        self.assertEqual(obs,exp_html)
        self.assertTrue(exists(filename1))
        self.assertTrue(exists(filename3))
        self.assertTrue(exists(filename4))
        self.assertTrue(exists(folder1)) 
        self.assertTrue(exists(folder2)) 
        self.assertTrue(exists(folder3)) 
        
    def test_save_rarefaction_data(self):
        '''save_rarefaction_data: This tests the rarefaction average output'''
        
        exp=['# test.txt\n', '# SampleID\n', 'xaxis: 10.0\t\n', 'xmax: 140\n', '>> Sample1\n', 'color #7fff00\n', 'series ', '2.03172\t9.441785\t12.508435\t\n', 'error ', 'nan\tnan\tnan\t\n']
        
        obs=save_rarefaction_data(self.ave_seqs_per_sample1,self.data['xaxis'],\
                                  self.xmax,self.mapping_category,self.colors2,\
                                  'test.txt',self.data_colors,self.groups,\
                                  self.std_type)
                    
        self.assertEqual(obs,exp)
        
    def test_save_single_rarefaction_plots(self):
        '''save_single_rarefaction_plots: this generates a plot with raw
           rarefaction data on a plot and tests whether a file is generated'''
           
        filename1='/tmp/testcol_0_row_0_raw.png'
        self._paths_to_clean_up = [filename1]

        exp={'test': {'groups': {'SampleID': {'Sample1': {'groupcolor': '#ff0000', 'raw_link': 'html_plots/testcol_0_row_0_raw.png', 'groupsamples': ['Sample1']}}}, 'samples': {'Sample1': {'color': '#ff0000', 'link': 'html_plots/testcol_0_row_0.png'}}}}
                     
        obs=save_single_rarefaction_plots(self.sample_dict, \
                    self.imagetype,self.metric_name, \
                    self.data_colors, self.colors, '/tmp/', \
                    self.background_color, self.label_color, self.resolution, \
                    self.ymax,self.xmax, \
                    self.rarefaction_legend_mat,self.groups, \
                    self.mapping_category,'Sample1',self.mapping_lookup)
 
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
           
        groups={'M': ['123','234','345']}
        
        test = make_error_series(self.ave_seqs_per_sample, groups, \
                                 self.std_type)
        
        exp2={'M': [0.81346064945802188, 4.2928534982243063, \
                    5.4903030649080646]}
        
        self.assertEqual(test[0], self.exp_err_series_ave)
        self.assertEqual(test[1], exp2)

    def test_make_error_series_stderr(self):
        '''make_error_series_stderr: this tests whether the standard errors 
           were correctly calculated'''
        
        groups={'M': ['123','234','345']}

        test2 = make_error_series(self.ave_seqs_per_sample, groups, \
                                 'stderr')
                                 
        exp2={'M': [0.57520354146018038, 3.0355058192348001, \
                    3.8822305279657781]}
                    
        self.assertEqual(test2[0], self.exp_err_series_ave)
        self.assertEqual(test2[1], exp2)

    def test_make_plots(self):
        """make_plots: tests whether the average plots are generated and if
           dictionary for the html generation is properly formatted"""

        filename1='/tmp/test/testcol_0_row_0_ave.png'
        filename2='/tmp/test/testcol_0_row_0_raw.png'
        folder1='/tmp/test/'
        
        self._paths_to_clean_up = [filename1,filename2]
        self._folders_to_cleanup=[folder1]

        exp1={'SampleID': {'Sample1': {'test': {'ave': ['     7.000', '     2.052'], 'err': ['       nan', '     0.000']}}}}
        exp2={'test': {'groups': {'SampleID': {'Sample1': {'groupcolor': '#ff0000', 'raw_link': 'html_plots/testcol_0_row_0_raw.png', 'groupsamples': ['Sample1'], 'ave_link': 'html_plots/testcol_0_row_0_ave.png'}}}, 'samples': {'Sample1': {'color': '#ff0000', 'link': 'html_plots/testcol_0_row_0.png'}}}}
        
        create_dir('/tmp/test/',False)
        
        obs1,obs2 = make_plots(self.background_color,self.label_color, \
                          self.rare_data,self.ymax, self.xmax,'/tmp/test/', \
                          self.resolution, self.imagetype,self.groups,\
                          self.colors,self.data_colors,self.metric_name,\
                          self.labelname,self.rarefaction_data_mat, \
                          self.rarefaction_legend_mat,self.sample_dict, \
                          self.data_colors,self.colors2,self.mapping_lookup)
        
        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,exp2)
        self.assertTrue(exists(filename1))
        self.assertTrue(exists(filename2))
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
'''\n<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n<html>\n<head>\n  <meta http-equiv="content-type" content="text/html;">\n  <title>Rarefaction Curves</title>\n<style type="text/css">\ntd.data{font-size:10px;border-spacing:0px 10px;text-align:center;}\ntd.headers{font-size:12px;font-weight:bold;text-align:center;}\ntable{border-spacing:0px;}\n.removed{display:none;}\n.expands{cursor:pointer; cursor:hand;}\n.child1 td:first-child{padding-left: 3px;}\n</style>\n<script language="javascript" type="text/javascript">\n\nfunction show_hide_category(checkobject){\n    var imagetype=document.getElementById(\'imagetype\').value;\n    img=document.getElementById(checkobject.name.replace(\'_raw\'+imagetype,\'_ave\'+imagetype))\n    if (checkobject.checked==false){\n        img.style.display=\'none\';\n    }else{\n        img.style.display=\'\';\n    }\n}\n\nfunction reset_tree(){\n    var category=document.getElementById(\'category\').value;\n    var metric=document.getElementById(\'metric\').value;\n    var old_all_categories=document.getElementById(\'all_categories\');\n    var imagetype=document.getElementById(\'imagetype\').value;\n    cat_list=old_all_categories.value.split(\'$#!\')\n    if (metric!=\'\' && category != \'\'){\n    for (var i=1, il=cat_list.length; i<il; i++){\n        group=metric+category+cat_list[i]\n        main_class=metric+category\n        var exp_item=document.getElementById(group);\n        if (exp_item!=null){\n            if (exp_item.innerHTML==\'\\u25BC\'){\n                exp_item.innerHTML=\'\\u25B6\'\n                var rows=document.getElementsByName(group);\n                for (var j=0, jl=rows.length; j<jl; j++){\n                    rows[j].style.display="none";\n                }\n            }\n            var rows=document.getElementsByName(group+\'_raw\'+imagetype);\n            for (var j=0, jl=rows.length; j<jl; j++){\n                if (rows[j].checked==false){\n                    rows[j].checked=true;\n                }\n            }\n        }\n    }\n}\n}\n\nfunction changeMetric(SelObject){\n    var category=document.getElementById(\'category\');\n    var old_metric=document.getElementById(\'metric\');\n    var imagetype=document.getElementById(\'imagetype\').value;\n    var legend=document.getElementById(\'legend\');\n    var array=document.getElementById(\'all_categories\').value.split(\'$#!\')\n    var plots=document.getElementById(\'plots\');\n    plots.style.display=\'none\'\n    reset_tree();\n    if (category.value != \'\'){\n        legend.style.display="";\n        cat=SelObject.value+category.value\n        data_display=document.getElementsByName(cat)\n        for (var i=0, il=data_display.length; i<il; i++){\n            data_display[i].style.display="";\n        }\n        cat=old_metric.value+category.value\n        data_hide=document.getElementsByName(cat)\n        for (var i=0, il=data_hide.length; i<il; i++){\n            data_hide[i].style.display="none";\n        }\n        data_display=document.getElementsByName(category.value)\n        for (var i=0, il=data_display.length; i<il; i++){\n            data_display[i].style.display="";\n        }\n        new_cat=SelObject.value+category.value\n        plots.innerHTML=\'\'\n        for (var i=1, il=array.length; i<il; i++){\n            img=document.createElement(\'img\')\n            img.setAttribute(\'width\',"600px") \n            img.setAttribute(\'id\',array[i]+\'_ave\'+imagetype)\n            img.setAttribute(\'style\',\'position:absolute;z-index:0\')\n            img.setAttribute(\'src\',"./html_plots/"+SelObject.value+array[i]+\'_ave\'+imagetype)\n            plots.appendChild(img)\n        }\n        plots.style.display=\'\'\n    }\n    \nold_metric.value=SelObject.value;\n}\n\nfunction changeCategory(SelObject){\n    var old_category=document.getElementById(\'category\');\n    var metric=document.getElementById(\'metric\').value;\n    var imagetype=document.getElementById(\'imagetype\').value;\n    var legend=document.getElementById(\'legend\');\n    var plots=document.getElementById(\'plots\');\n    var array=SelObject.value.split(\'$#!\')\n    var old_all_categories=document.getElementById(\'all_categories\')\n    category=array[0]\n    plots.style.display=\'none\'\n    reset_tree();\n\n    if (metric != \'\'){\n        legend.style.display="";\n\n        data_display=document.getElementsByName(category)\n        for (var i=0, il=data_display.length; i<il; i++){\n            data_display[i].style.display="";\n        }\n        data_hide=document.getElementsByName(old_category.value)\n        for (var i=0, il=data_hide.length; i<il; i++){\n            data_hide[i].style.display="none";\n        }\n        cat=metric+category\n        data_display=document.getElementsByName(cat)\n        for (var i=0, il=data_display.length; i<il; i++){\n            data_display[i].style.display="";\n        }\n        cat=metric+old_category.value\n        data_hide=document.getElementsByName(cat)\n        for (var i=0, il=data_hide.length; i<il; i++){\n            data_hide[i].style.display="none";\n        }\n        cat=metric+category\n        plots.innerHTML=\'\'\n        for (var i=1, il=array.length; i<il; i++){\n            img=document.createElement(\'img\')\n            img.setAttribute(\'width\',"600px") \n            img.setAttribute(\'id\',metric+array[i]+\'_ave\'+imagetype)\n            img.setAttribute(\'style\',\'position:absolute;z-index:0\')\n            img.setAttribute(\'src\',"./html_plots/"+metric+array[i]+\'_ave\'+imagetype)\n            plots.appendChild(img)\n        }\n        plots.style.display=\'\'\n    }\nold_all_categories.value=SelObject.value;\nold_category.value=category;\n}\nfunction toggle(){\n    var plots=document.getElementById(\'plots\');\n    var imagetype=document.getElementById(\'imagetype\').value;\n    var plot_str=\'\';\n    var category=document.getElementById(\'category\');\n    var metric=document.getElementById(\'metric\');\n    expansion_element=document.getElementById(arguments[0]);\n    rows=document.getElementsByName(arguments[0]);\n    if (expansion_element.innerHTML==\'\\u25B6\'){\n        expansion_element.innerHTML=\'\\u25BC\'\n        show_row=arguments[0]+\'_raw\'+imagetype\n        \n        if (document.getElementById(show_row)==null){\n            img=document.createElement(\'img\')\n            img.setAttribute(\'width\',"600px") \n            img.setAttribute(\'id\',arguments[0]+\'_raw\'+imagetype)\n            img.setAttribute(\'style\',\'position:absolute;z-index:0\')\n            img.setAttribute(\'src\',"./html_plots/"+arguments[0]+\'_raw\'+imagetype)\n            plots.appendChild(img)\n        }else{\n            document.getElementById(arguments[0]+\'_raw\'+imagetype).style.display=\'\'\n        }\n        for (var i=0, il=rows.length;i<il;i++){\n            rows[i].style.display=\'\';\n        }\n    }else{\n        expansion_element.innerHTML=\'\\u25B6\'\n        document.getElementById(arguments[0]+\'_raw\'+imagetype).style.display=\'none\'\n        for (var i=0, il=rows.length;i<il;i++){\n            rows[i].style.display=\'none\';\n        }\n    }\n}\n\nfunction show_hide_categories(SelObject){\n    var all_categories=document.getElementById(\'all_categories\').value.split(\'$#!\')\n    var category=document.getElementById(\'category\').value;\n    var imagetype=document.getElementById(\'imagetype\').value;\n    var metric=document.getElementById(\'metric\').value;\n    for (var i=1, il=all_categories.length; i<il; i++){\n        basename=metric+category+all_categories[i]\n        raw_image=basename+\'_raw\'+imagetype\n        ave_image=basename+\'_ave\'+imagetype\n        checkbox=document.getElementsByName(raw_image)\n        if (SelObject.value==\'All\'){\n            if (checkbox[0].checked==false){\n                checkbox[0].checked=true\n                document.getElementById(ave_image).style.display=\'\'\n            }\n        }else if (SelObject.value==\'None\'){\n            if (checkbox[0].checked==true){\n                checkbox[0].checked=false\n                document.getElementById(ave_image).style.display=\'none\'\n            }\n        }else if (SelObject.value==\'Invert\'){\n            if (checkbox[0].checked==true){\n                checkbox[0].checked=false\n                document.getElementById(ave_image).style.display=\'none\'\n            }else if (checkbox[0].checked==false){\n                checkbox[0].checked=true\n                document.getElementById(ave_image).style.display=\'\'\n            }\n        }\n    }\n    document.getElementById(\'show_category\').selectedIndex=0;\n}\n</script>\n\n</head>\n<body>\n<form action=\'\'>\n<input id="metric" type="hidden">\n<input id="category" type="hidden">\n<input id="imagetype" type="hidden" value=".png">\n<input id="all_categories" type="hidden">\n</form>\n<table><tr>\n<td><b>Select a Metric:</b></td>\n<td>\n<select onchange="javascript:changeMetric(this)">\n<option>&nbsp;</option>\n<option value="test">test</option>\n</select>\n</td>\n<td><b>&nbsp;&nbsp;Select a Category:</b></td>\n<td>\n<select onchange="javascript:changeCategory(this)">\n<option>&nbsp;</option>\n<option value="SampleID$#!col_0_row_0">SampleID</option>\n</select>\n</td>\n</table>\n<br>\n<div style="width:950px">\n<div id="plots" style="width:650px;height:550px;float:left;"></div>\n\n<div id="legend" style="width:300px;height:550px;float:right;display:none;">\n    <p><b>Show Categories: \n    <select id="show_category" onchange="show_hide_categories(this);">\n        <option value="">&nbsp;</option>\n        <option value="All">All</option>\n        <option value="None">None</option>\n        <option value="Invert">Invert</option>\n    </select>\n    </b></p>\n<b>Legend</b><div STYLE="border: thin black solid; height: 300px; width: 200px; font-size: 12px; overflow: auto;"><table>\n<tr id="testSampleID" name="testSampleID" style="display: none;"><td class="data" onmouseover="document.body.style.cursor=\'pointer\'"  onmouseout="document.body.style.cursor=\'default\'" onclick="toggle(\'testcol_0_row_0\')" id="testcol_0_row_0" name="\'Sample1\'">&#x25B6;</td><td><input name="testcol_0_row_0_raw.png" type="checkbox" checked="True" onclick="show_hide_category(this)"></td><td style="color:#ff0000">&#x25A0;&nbsp;</td><td class="data"><b>Sample1</b></td></tr>\n<tr id="testcol_0_row_0_raw" name="testcol_0_row_0" style="display: none;"><td class="data" align="right">&#x221F;</td><td></td><td style="color:#ff0000">&#x25C6;</td><td class="data" align="left"><b>Sample1</b></td></tr>\n</table></div></div>\n<div style="position:relative;clear:both;">\n<table id="rare_data" border="1px">\n<tr name="SampleID" style="display: none;"><td class="headers">SampleID</td><td class="headers">Seqs/Sample</td>\n<td class="headers">test Ave.</td><td class="headers">test Err.</td>\n</tr>\n<tr name="SampleID" style="display: none;">\n<td class="data" bgcolor="#ff0000">Sample1</td><td class="data">10.0</td>\n<td class="data">     7.000</td><td class="data">       nan</td>\n</tr>\n</table>\n</div>\n</div>\n</body>\n</html>\n'''


if __name__ == "__main__":
    main()
