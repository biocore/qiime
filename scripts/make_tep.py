#!/usr/bin/env python
# File created on 15 Jun 2011
from __future__ import division

__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Meg Pirrung"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Meg Pirrung"
__email__ = "meganap@gmail.com"
__status__ = "Development"

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
 make_option('-p','--prefs_file_dir',type="string",help='Path to prefs file', dest='prefs_file_fp'),\

 # make_option('--create_jnlp',action='store_true',
 #  help='create a jnlp file [default: %default]'),\
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
	
	<jar href="topiaryexplorer0.9.4.jar" />
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

data_color_hsv = {
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

def make_te_prefs(prefs_dict):
    sample_coloring = prefs_dict['sample_coloring']
    lines = []
    for k in sample_coloring:
        for t in sample_coloring[k]['colors']:
            if(type(t) == tuple):
                lines.append(''.join([str(i)+',' for i in t[1]])+'\n')
            if(type(t) == str):
                lines.append(t+':'+''.join([str(i)+',' for i in data_color_hsv[sample_coloring[k]['colors'][t]]])+'\n')
        lines.append('>default'+k+':'+k+'\n')
    # print lines
    return lines

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
    
    if opts.prefs_file_fp:
        prefs_fp = opts.prefs_file_fp
        prefs_dict = eval(open(prefs_fp,'U').read())
        te_prefs = make_te_prefs(prefs_dict)
        lines += '\n>>pre\n'
        lines += te_prefs
    
    tepfile.writelines(lines)

    # if opts.create_jnlp:
    jnlpfile = open(jnlp_fp, 'w')
    lines = [jnlp_top_block]
    if(opts.web_flag):
        lines += 'http://topiaryexplorer.sourceforge.net/app/'
    else:
        topiaryexplorer_project_dir =\
         load_qiime_config()['topiaryexplorer_project_dir']
        if topiaryexplorer_project_dir == None:
            option_parser.error("Couldn't create jnlp file - topiaryexplorer_project_dir is not defined in your qiime_config. tep file was created sucessfully.")
        lines += 'file:' + topiaryexplorer_project_dir
    
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