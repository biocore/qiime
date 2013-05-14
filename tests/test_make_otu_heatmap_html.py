#!/usr/bin/env python

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project" #consider project name
__credits__ = ["Jesse Stombaugh", "Jose Carlos Clemente Litran",
               "Jai Ram Rideout"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"

from numpy import array, log
import shutil
from shutil import rmtree
from os.path import join
from cogent.util.unit_test import TestCase, main
from qiime.make_otu_heatmap_html import (make_html_doc,create_javascript_array,
                                         filter_by_otu_hits,
                                         get_log_transform,
                                         generate_heatmap_plots)
from qiime.util import create_dir,get_qiime_project_dir
from biom.table import table_factory

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        self.output_dir='/tmp/'
        
        otu_table_vals = array([[0,0],[1,5]])

        self.otu_table = table_factory(otu_table_vals,
                                        ['Sample1', 'Sample2'],
                                        ['OTU1', 'OTU2'],
                                        [None, None],
                                        [{"taxonomy": ["Bacteria"]},
                                         {"taxonomy": ["Archaea"]}])

        filt_otu_table_vals = array([[1,5]])

        self.filt_otu_table = table_factory(filt_otu_table_vals,
                                             ['Sample1', 'Sample2'],
                                             ['OTU2'],
                                             [None, None],
                                             [{"taxonomy": ["Archaea"]}])

        self.num_otu_hits=5
        self._folders_to_cleanup=[]
        
    def tearDown(self):
        '''This function removes the generated files'''
        map(rmtree,self._folders_to_cleanup)
        
    def test_make_html_doc(self):
        """make_html_doc: create the web interface for otu heatmaps"""
        obs = make_html_doc('./test/test_otu_count_table.js')

        self.assertEqual(obs,exp_html_script)
            
    def test_create_javascript_array(self):
        """create_javascript_array: takes a numpy array and generates a
javascript array"""
        obs = create_javascript_array(self.filt_otu_table)

        self.assertEqual(obs,exp_js_array)        
            
    def test_filter_by_otu_hits(self):
        """filter_by_otu_hits: filters the table by otu hits per otu"""
        exp = array([['#OTU ID', 'OTU2'],['Sample1',1],['Sample2',5],\
                         ['Consensus Lineage','Archaea;']])
        
        obs = filter_by_otu_hits(self.num_otu_hits, self.otu_table)

        # see note in test_get_log_transform about this assert
        self.assertEqual(obs._data.items(), self.filt_otu_table._data.items())
    
    def test_get_log_transform(self):
        orig_data = array([[0,1,2],[1000,0,0]])

        orig_otu_table = table_factory(orig_data,
                                        ['Sample1', 'Sample2', 'Sample3'],
                                        ['OTU1', 'OTU2'],
                                        [None, None, None],
                                        [{"taxonomy": ["Bacteria"]},
                                         {"taxonomy": ["Archaea"]}])

        exp_data = array([[0,0.69314718,1.38629436],[7.60090246,0,0]])
        exp_otu_table = table_factory(exp_data,
                                       ['Sample1', 'Sample2', 'Sample3'],
                                       ['OTU1', 'OTU2'],
                                       [None, None, None],
                                       [{"taxonomy": ["Bacteria"]},
                                        {"taxonomy": ["Archaea"]}])

        log_otu_table = get_log_transform(orig_otu_table, eps = None)

        # comparing directly log_otu_table against exp_otu_table doesn't work,
        #  needs to be modified in the otu table object
        self.assertFloatEqual(list(log_otu_table.iterSampleData()),
                              list(exp_otu_table.iterSampleData()))

    def test_generate_heatmap_plots(self):
        """generate_heatmap_plots: create default output files"""
        
        # create directories and move js files to verify everything works
        # in the script file
        dir_path=join(self.output_dir,'test')
        create_dir(dir_path)
        
        js_dir_path=join(dir_path,'js')
        create_dir(js_dir_path)
        
        self._folders_to_cleanup.append(dir_path)
        
        qiime_dir=get_qiime_project_dir()
        
        js_path=join(qiime_dir,'qiime/support_files/js')
        shutil.copyfile(join(js_path,'overlib.js'), join(js_dir_path,'overlib.js'))
        shutil.copyfile(join(js_path,'otu_count_display.js'), join(js_dir_path,'otu_count_display.js'))
        shutil.copyfile(join(js_path,'jquery.js'), join(js_dir_path,'jquery.js'))
        shutil.copyfile(join(js_path,'jquery.tablednd_0_5.js'), join(js_dir_path,'jquery.tablednd_0_5.js'))
        
    
        # generate otu_table object
        orig_data = array([[0,1,2],[1000,0,0]])
                     
        orig_otu_table = table_factory(orig_data,
                                        ['Sample1', 'Sample2', 'Sample3'],
                                        ['OTU1', 'OTU2'],
                                        [None, None, None],
                                        [{"taxonomy": ["Bacteria"]},
                                         {"taxonomy": ["Archaea"]}])
        
        # put in an OTU sort order and sample order
        otu_sort=['OTU2','OTU1']
        sample_sort=['Sample2','Sample1','Sample3']
        num_otu_hits=3
        
        # generate test files
        generate_heatmap_plots(num_otu_hits, orig_otu_table, otu_sort, 
                               sample_sort, dir_path, js_dir_path, 
                               'test',fractional_values=False)
        
        
        self.assertEqual(open(join(js_dir_path,'test.js'),'U').read(), 
                         exp_js_output_file)
        
exp_js_array='''\
var OTU_table=new Array();\n\
var i=0;\n\
for (i==0;i<4;i++) {\n\
OTU_table[i]=new Array();}\n\
OTU_table[0][0]='#OTU ID';\n\
OTU_table[0][1]='OTU2';\n\
OTU_table[1][0]='Sample1';\n\
OTU_table[1][1]=1;\n\
OTU_table[2][0]='Sample2';\n\
OTU_table[2][1]=5;\n\
OTU_table[3][0]='Consensus Lineage';\n\
OTU_table[3][1]='Archaea';\n\
'''

exp_js_output_file="""\
var otu_num_cutoff=3;
var OTU_table=new Array();
var i=0;
for (i==0;i<5;i++) {
OTU_table[i]=new Array();}
OTU_table[0][0]='#OTU ID';
OTU_table[0][1]='OTU2';
OTU_table[0][2]='OTU1';
OTU_table[1][0]='Sample2';
OTU_table[1][1]=0;
OTU_table[1][2]=1;
OTU_table[2][0]='Sample1';
OTU_table[2][1]=1000;
OTU_table[2][2]=0;
OTU_table[3][0]='Sample3';
OTU_table[3][1]=0;
OTU_table[3][2]=2;
OTU_table[4][0]='Consensus Lineage';
OTU_table[4][1]='Archaea';
OTU_table[4][2]='Bacteria';
"""

exp_html_script = \
    r'''
    <html>
    <head>
    	<script type="text/javascript" src="js/overlib.js"></script>
        <script type="text/javascript" src="./test/test_otu_count_table.js"></script>
    	<script type="text/javascript" src="js/otu_count_display.js"></script>
    	<script type="text/javascript" src="./js/jquery.js"></script>
    	<script type="text/javascript" src="./js/jquery.tablednd_0_5.js"></script>
        <script type="text/javascript">

    
        $(document).ready(function() {
    
        	$('#otu_table_body').tableDnD({
        		onDragStart: function(table, new_row) {
        			if (row==new_row.parentNode.rowIndex && is_selected==1){
        				change_sel_row=1;
        			}else{
        				old_row=new_row.parentNode.rowIndex;
        				change_sel_row=0;
        			}
        		},
        		onDrop: function(table, new_row) {
        			if (change_sel_row==1){
        				row=new_row.rowIndex;
        			}else if(old_row<row && new_row.rowIndex>row){
        				row=row-1;
        			}else if(old_row>row && new_row.rowIndex<row){
        				row=row+1;
        			}
        		},
            	dragHandle: "dragHandle"
        	});
            var otu_cutoff=document.getElementById("otu_count_cutoff");
            otu_cutoff.value=otu_num_cutoff;    
        });
        </script>
    	<style type="text/css">
    	    th.rotate{ 
    			white-space : nowrap;
    			-webkit-transform: rotate(-90deg) translate(20px, 0px); 
    			-moz-transform: rotate(-90deg) translate(20px, 0px);	
    			font-family:arial;
    			font-size:9px;
    		}
    		th.lineage{ 
        	    white-space : nowrap;
        	    text-align:left;
        	    font-family:arial;
        	    font-size:10px;
        	    font-weight: bolder;
        	}
        	td.dragHandle{ 
            	white-space : nowrap;
            	text-align:left;
            	font-family:arial;
            	font-size:10px;
            	font-weight: bolder;
        	}
        	td{ 
            	white-space : nowrap;
            	font-family:arial;
            	font-size:10px;
            	text-align:center;
            	font-weight: bolder;
        	}       
        	table{ 
            	border-spacing: 0;
            	text-align:center;
        	}
        	p{
            		text-align:left;
            		font-weight: normal;
        	}    
    	</style>
    </head>
    <body>
    	<p>
    		Filter by Counts per OTU: <input type="text" id="otu_count_cutoff" value="">
    		<input type="button" onclick="javascript:create_OTU_intervals();" value="Sample ID">
    		<input type="button" onclick="javascript:write_taxon_heatmap();" value="Taxonomy">
    	</p>
    	<br><br><br><br><br><br>
    	<table id='otu_table_html'>
    		<thead id='otu_table_head'>
    		</thead>
    		<tbody id='otu_table_body'>
    		<tr><td class="dragHandle"></td>
    		</tr>
    		<tr><td class="dragHandle"></td>
    		</tr>
    		</tbody>
    	</table>

    </body>
    </html>'''

#run tests if called from command line
if __name__ == "__main__":
    main()
