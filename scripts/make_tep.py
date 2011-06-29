#!/usr/bin/env python
# File created on 15 Jun 2011
from __future__ import division

__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Meg Pirrung"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Meg Pirrung"
__email__ = "meganap@gmail.com"
__status__ = "Release"

from os.path import split,splitext
from qiime.util import parse_command_line_parameters, make_option
from qiime.util import load_qiime_config, create_dir
from qiime.parse import parse_otu_table
import os

script_info = {}
script_info['brief_description'] = "Makes TopiaryExplorer project file"
script_info['script_description'] = "This script makes a TopiaryExplorer project file (.tep) and a jnlp file with the data location preloaded.\n\nWARNING: The jnlp file relies on an absolute path, if you move the .tep file, the generated jnlp will no longer work. However, you can still open the .tep file from your normal TopiaryExplorer install."
script_info['script_usage'] = []
script_info['script_usage'].append(("Example","Create .tep file and .jnlp file:","make_tep.py -i otu_table.txt -m mapping_file.txt -t otus.tre -o my_data"))
# script_info['script_usage'].append(("""""","Create .tep file and .jnlp file, then launch TopiaryExplorer:","make_tep.py -i otu_table.txt -m mapping_file.txt -t otus.tre -o my_data -l"))
script_info['output_description']= "The result of this script is written to a .tep file and a .jnlp file, both with the name supplied by -o"
script_info['required_options'] = [\
 # Example required option
 make_option('-i','--otu_file',type="string",help='Path to read otu table',dest='otu_table_fp'),\
 make_option('-m','--mapping_file',type="string",help='Path to read mapping file',dest='mapping_fp'),\
 make_option('-t','--tree_file',type="string",help='Path to read tree',dest='tree_fp')
]
script_info['optional_options'] = [\
 # Example optional option
 make_option('-o','--output_dir',type="new_dirpath",help='the output directory [default: %default]', dest='out_fp'),\
  make_option('-w','--web',action='store_true',default=False, help='web codebase jnlp flag [default: %default]', dest='web_flag'),
  make_option('-u','--url_path',type="string",help='url path', dest='url')
 # make_option('-l','--launch_TE',action='store_false', default="false",help='Option to launch TopiaryExplorer [default: %default]', dest='launch')
]
script_info['version'] = __version__

jnlp_top_block = """<?xml version="1.0" encoding="utf-8"?>

<jnlp codebase=\""""

jnlp_middle_block = """\">
    
    <information>
    <title>TopiaryExplorer</title>
    <vendor>University of Colorado</vendor>
    <description>TopiaryExplorer</description>
    
    <offline-allowed/>
     

  </information>
  
  <security>
      <all-permissions/>
  </security>
  
  <resources>
    <j2se version="1.6+" initial-heap-size="500M" max-heap-size="2000m" />
	<extension name="jogl" href="http://download.java.net/media/jogl/builds/archive/jsr-231-webstart-current/jogl.jnlp" />
	
	<jar href="topiaryexplorer.jar" />
	<jar href="lib/core.jar" />
	<jar href="lib/itext.jar" />
	<jar href="lib/pdf.jar" />
	<jar href="lib/ojdbc14.jar" />
	<jar href="lib/opengl.jar" />
	<jar href="lib/mysql-connector-java-5.1.10-bin.jar" />
	<jar href="lib/javaws.jar" />
	<jar href="lib/classes12.jar" />
	<jar href="lib/jogl.jar" />
	
  </resources>
  
  <application-desc main-class="topiaryexplorer.TopiaryExplorer">
  <argument>"""

jnlp_bottom_block = """</argument>
	</application-desc>
</jnlp>
"""

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    otu_table_fp = opts.otu_table_fp
    mapping_fp = opts.mapping_fp
    tree_fp = opts.tree_fp
    output_dir = opts.out_fp
    output_basename = splitext(split(otu_table_fp)[1])[0]
    
    if not output_dir:
        output_dir = 'make_tep_output/'
    
    create_dir(output_dir)
    
    tep_fp = '%s/%s.tep' % (output_dir,output_basename)      # opts.out_fp+'.tep'
    jnlp_fp = '%s/%s.jnlp' % (output_dir,output_basename)
    tepfile = open(tep_fp, 'w')
    otu_lines = open(otu_table_fp, 'U').readlines()
    sample_ids, otu_ids, otu_table, metadata = parse_otu_table(otu_lines)
    mapping_lines = open(mapping_fp, 'U')    
    tree_lines = open(tree_fp, 'U')
    
    lines = ['>>tre\n']
    lines += tree_lines.readlines() 
    lines += '\n'
    if(metadata):
        lines += '>>otm\n#OTU ID\tOTU Metadata\n'
        for i in range(len(otu_ids)):
            lines += otu_ids[i] + '\t'
            for m in metadata[i]:
                lines += m + ';'
            # lines = lines[:len(lines)-1]
            lines += '\n'
    lines += '>>osm\n'
    lines += otu_lines
    lines += '\n>>sam\n'
    lines += mapping_lines.readlines()
    
    tepfile.writelines(lines)
    
    jnlpfile = open(jnlp_fp, 'w')
    lines = [jnlp_top_block]
    if(opts.web_flag):
        lines += 'http://topiaryexplorer.sourceforge.net/app/'
    else:
        lines += 'file:'+load_qiime_config()['topiaryexplorer_project_dir']
    lines += jnlp_middle_block
    if(opts.url):
        lines += opts.url
    else:
        lines += os.path.abspath(tep_fp)
    # lines += os.path.abspath(tep_fp)
    lines += jnlp_bottom_block
    jnlpfile.writelines(lines)
    
if __name__ == "__main__":
    main()