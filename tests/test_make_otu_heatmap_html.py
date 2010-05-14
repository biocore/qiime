#!/usr/bin/env python
#file test_make_otu_heatmap_html.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"

from numpy import array
from cogent.util.unit_test import TestCase, main
from qiime.make_otu_heatmap_html import (make_html_doc,create_javascript_array,\
                                       filter_by_otu_hits,line_converter)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        self.col_header=['Sample1', 'Sample2']
        self.row_header=['OTU1','OTU2']
        self.otu_table=array([[0,0],[1,5]])
        self.lineages=[['Bacteria'],['Archaea']]

        self.data={}
        self.data['otu_counts']=self.col_header,self.row_header,self.otu_table,\
                                self.lineages

        self.num_otu_hits=5
            
    def test_make_html_doc(self):
        """make_html_doc: create the web interface for otu heatmaps"""
        obs = make_html_doc('./test/test_otu_count_table.js')

        self.assertEqual(obs,exp_html_script)
            
    def test_create_javascript_array(self):
        """create_javascript_array: takes a numpy array and generates a
javascript array"""
        self.rows=[['#OTU ID', 'OTU2'],['Sample1',1],['Sample2',5],\
                   ['Consensus Lineage','Archaea']]

        obs=create_javascript_array(self.rows)

        self.assertEqual(obs,exp_js_array)        
            
    def test_filter_by_otu_hits(self):
        """filter_by_otu_hits: filters the table by otu hits per otu"""
        exp=array([['#OTU ID', 'OTU2'],['Sample1',1],['Sample2',5],\
             ['Consensus Lineage','Archaea;']])
        
        obs=filter_by_otu_hits(self.num_otu_hits,self.data)
 
        self.assertEqual(obs,exp)

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
