#!/usr/bin/env python
#file make_otu_heatmap_html.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"


from numpy import array,concatenate,asarray,transpose
from cogent.parse.table import SeparatorFormatParser
from optparse import OptionParser
from qiime.parse import parse_otu_table
import os

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
                js_array+="OTU_table[%i][%i]='%s';\n" % (i,j,(rows[i][j]))
            else:
                js_array+="OTU_table[%i][%i]=%s;\n" % (i,j,int(rows[i][j]))
            
    return js_array

def filter_by_otu_hits(num_otu_hits,data):
    """Filter the OTU table by the number of otus per sample"""
    col_header,row_header,otu_table,lineages=data['otu_counts']
    lineages_update=[]
    for i in range(len(lineages)):
        new_lineages=''
        for j in lineages[i]:
            new_lineages+=j.strip('"')+';'
        lineages_update.append(new_lineages)
    lineages_update=array(lineages_update)
    rows_filtered=[]
    row_head=concatenate((['#OTU ID'],col_header))
    lineage_head=array(['Consensus Lineage'])
    row_head2=concatenate((row_head,lineage_head))
    rows_filtered.append(row_head2)

    for i in range(len(otu_table)):
        if otu_table[i].sum()>num_otu_hits:
            row_head_otu_count=concatenate(([row_header[i]],otu_table[i],[lineages_update[i]]))

            rows_filtered.append(row_head_otu_count)

    rows_filtered=array(rows_filtered) 
    trans_rows_filtered=rows_filtered.transpose()
    return trans_rows_filtered
    
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
    
def get_otu_counts(fpath, data):
    """Reads the OTU table file into memory"""
    try:
        sample_ids,otu_ids,otu_table,lineages=parse_otu_table(open(fpath))    
    except (TypeError, IOError):
        raise MissingFileError, 'OTU Count file required for this analysis'
    
    if lineages==[]:
        raise ValueError, '\n\nThe lineages are missing from the OTU table.  If you used single_rarefaction.py to create your otu_table, make sure you pass the "--lineages_included" option.\n'
        
    return sample_ids,otu_ids,otu_table,lineages

def generate_heatmap_plots(options,data, dir_path, js_dir_path,filename):
    """Generate HTML heatmap and javascript array for OTU counts"""

    #Convert number of otu hits argument into an integer
    num_otu_hits=int(options.num_otu_hits)
    #Filter by number of OTU hits
    rows=filter_by_otu_hits(num_otu_hits, data)
    #print rows
    #print rows[0]
    
    # This sorts the otus by the tree supplied
    if data['otu_order']:
        #print data['otu_order']
        new_otu_table=[]
        tran_rows=rows.transpose()
    
        new_otu_table.append(tran_rows[0])
        for i in data['otu_order']:
            for j in tran_rows:
                if i==j[0]:
                    new_otu_table.append(j)
        rows= asarray(new_otu_table).transpose()
        
    #Convert OTU counts into a javascript array
    js_array=create_javascript_array(rows)

    #Write otu filter number
    js_otu_cutoff='var otu_num_cutoff=%i;' % num_otu_hits
    
    #Write js array to file
    js_filename=os.path.join(js_dir_path,filename)+'.js'
    jsfile = open(js_filename,'w')
    jsfile.write(js_otu_cutoff)
    jsfile.write(js_array)
    jsfile.close()

    #Write html file
    html_filename=os.path.join(dir_path,filename)+'.html'
    js_file_location='js/'+filename+'.js'
    table_html=make_html_doc(js_file_location)
    ofile = open(html_filename,'w')
    ofile.write(table_html)
    ofile.close()
