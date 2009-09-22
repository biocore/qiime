#!/usr/bin/env python
#file gen_otu_heatmap_html.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2009, the 454 Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Prototype"

"""
Author: Jesse Stombaugh (jesse.stombaugh@colorado.edu)
Status: Prototype

Requirements:
Python 2.5

Example 1: Create html file and javascript array from otu counts table:
Usage: python gen_otu_heatmap_html.py -o otu_counts.txt

Example 2: Create html file, then javascript array in specified directory:
Usage: python gen_otu_heatmap_html.py -o otu_counts.txt -x ./test/

Example 3: Create html file, then javascript array where the number of hits
per otu are speified:
Usage: python gen_otu_heatmap_html.py -o otu_counts.txt -n 50

"""

import numpy
from numpy import array
from cogent.parse.table import SeparatorFormatParser
from optparse import OptionParser
from gen_3d_plots import create_dir
import shutil

def make_html_doc(js_filename):
    """Create the basic framework for the OTU table heatmap"""
    html_script = \
    r'''
    <html>
    <head>
    	<script type="text/javascript" src="js/overlib.js"></script>
        <script type="text/javascript" src="%s"></script>
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
    </html>''' % (js_filename)
    return html_script

def create_javascript_array(rows):
    """Convert the OTU table counts into a javascript array"""
    
    js_array='\
    var OTU_table=new Array();\n\
    var i=0;\n\
    for (i==0;i<%i;i++) {\n\
    OTU_table[i]=new Array();}\n' % (len(rows))
    
    for i in range(len(rows)):
        for j in range(len(rows[i])):
            if i==0 or j==0 or i==len(rows)-1:
                js_array+='OTU_table[%i][%i]="%s";\n' % (i,j,(rows[i][j]))
            else:
                js_array+='OTU_table[%i][%i]=%i;\n' % (i,j,(rows[i][j]))
            
    return js_array

def filter_by_otu_hits(num_otu_hits,data):
    """Filter the OTU table by the number of otus per sample"""
    rows=data['otu_counts']

    rows_filtered=[]
    for i in range(len(rows)):
        if (rows[i,1:len(rows[i])-1].sum())>num_otu_hits:
            rows_filtered.append(rows[i])
    rows_filtered=array(rows_filtered)
    rows = rows_filtered.transpose()

    return rows
    
def line_converter():
    """Converts line elements into int's if possible"""
    def callable(line):
        new = []
        append = new.append
        for element in line:
            try:
                append(int(element))
            except ValueError:
                append(element)
        return new
    return callable
    
def get_otu_counts(options, data):
    """Reads the OTU table file into memory"""
    try:
        reader = SeparatorFormatParser(converter=line_converter(), sep='\t')
        rows = [line for line in reader(open(options.otu_count_fname))]
        title = rows.pop(0)[0]
    except (TypeError, IOError):
        raise MissingFileError, 'OTU Count file required for this analysis'

    data['otu_counts'] = numpy.asarray(rows, dtype='O')

    return data['otu_counts']

def _make_cmd_parser():
    """Returns the command-line options"""
    parser = OptionParser(usage="Usage: this_file.py -o <otu count output files>\
-n <num of otu hits: default=5> -x <write to directory: default=random dir>")
    parser.add_option('-o', '--otu_count', dest='otu_count_fname', \
        help='name of otu count file')
    parser.add_option('-x', '--dir-prefix', dest='dir_path',\
        help='directory prefix for all analyses')
    parser.add_option('-n', '--num_otu_hits', dest='num_otu_hits',\
        help='number of hits per OTU')
    options, args = parser.parse_args()
    return options

def _process_prefs(options):
    """opens files as necessary based on prefs"""
    data = {}
    
    if not options.num_otu_hits:
        options.num_otu_hits=5

    #Open and get coord data
    data['otu_counts'] = get_otu_counts(options, data)

    filepath=options.otu_count_fname
    filename=filepath.strip().split('/')[-1].split('.')[0]
    
    dir_path = create_dir(options.dir_path)
    js_dir_path = create_dir(dir_path+'js/')
    
    file_path=__file__.split('/')
    if len(file_path)==1:
       qiime_dir='./'
    else:
       qiime_dir='';
       for i in range(len(file_path)-1):
           qiime_dir+=file_path[i]+'/'
    
    shutil.copyfile(qiime_dir+'/js_library/overlib.js', js_dir_path+'overlib.js')
    shutil.copyfile(qiime_dir+'/js_library/otu_count_display.js', js_dir_path+\
                    'otu_count_display.js')
    shutil.copyfile(qiime_dir+'./js_library/jquery.js', js_dir_path+'jquery.js')
    shutil.copyfile(qiime_dir+'./js_library/jquery.tablednd_0_5.js', js_dir_path+\
                    'jquery.tablednd_0_5.js')

    
    
    action_str = '_do_heatmap_plots'
    try:
        action = eval(action_str)
    except NameError:
        action = None
    #Place this outside try/except so we don't mask NameError in action
    if action:
        action(options,data, dir_path,js_dir_path,filename)

def _do_heatmap_plots(options,data, dir_path, js_dir_path,filename):
    """Generate HTML heatmap and javascript array for OTU counts"""

    #Convert number of otu hits argument into an integer
    num_otu_hits=int(options.num_otu_hits)
    #Filter by number of OTU hits
    rows=filter_by_otu_hits(num_otu_hits, data)

    #Convert OTU counts into a javascript array
    js_array=create_javascript_array(rows)

    #Write otu filter number
    js_otu_cutoff='var otu_num_cutoff=%i;' % num_otu_hits
    
    #Write js array to file
    js_filename=js_dir_path+filename+'.js'
    jsfile = open(js_filename,'w')
    jsfile.write(js_otu_cutoff)
    jsfile.write(js_array)
    jsfile.close()

    #Write html file
    html_filename=dir_path+filename+'.html'
    js_file_location='js/'+filename+'.js'
    table_html=make_html_doc(js_file_location)
    ofile = open(html_filename,'w')
    ofile.write(table_html)
    ofile.close()

if __name__ == "__main__":
    from sys import argv, exit
    options = _make_cmd_parser()
    
    _process_prefs(options)
