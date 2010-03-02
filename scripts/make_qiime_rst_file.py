#!/usr/bin/env python
# File created on 15 Feb 2010
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Pre-release"


from qiime.util import parse_command_line_parameters
from optparse import make_option
from subprocess import Popen, PIPE, STDOUT
import os
from string import replace
import types
import re
from sys import exit, stderr

rst_text= \
'''\
.. _%s:

.. index:: %s

:mod:`%s` -- %s
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

%s


**Usage:** :file:`%s [options]`

**Input Arguments:**

.. note::

%s

**Output:**

%s

%s

'''


script_description = """This script will take a script file and convert the \
usage strings and options to generate a documentation .rst file."""

script_usage = """\
Example:
$ make_qiime_rst_file.py -i make_2d_plots.py -o doc/
"""

required_options = [\
 # Example required option
 make_option('-i','--input_script',help='This is the input script for which to \
 make a .rst file'),
 make_option('-o','--output_dir',help='the output directory [default: %default]')
]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]')
]



def convert_py_file_to_link(input_str):
    m=re.compile('[\w]+\.py')
    python_script_names=m.findall(input_str)
    
    if python_script_names:
        for i in python_script_names:
            individual_script_name=os.path.splitext(i)
            script_w_link=m.sub('`'+individual_script_name[0] + ' <./' + individual_script_name[0] + '.html>`_', input_str)
        
        return script_w_link
    else:
        return input_str



def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)


    #Determine if the input is a directory containing python scripts or a single
    #file and then create a list of those scripts
    if os.path.isdir(opts.input_script):
        file_names = os.listdir(opts.input_script)
        #Only take files which do not start with "." and end with a ".py"
        file_names = [fname for fname in file_names if not fname.startswith('.') and fname.endswith('.py')]
        
    elif os.path.isfile(opts.input_script):
        file_names=[str(opts.input_script)]
    
    else:
        print("io error, input path not valid.  Does it exist?")
        exit(1)

    
    #Identify the directory where results should be written.
    dir_path = opts.output_dir
    if dir_path and not dir_path.endswith('/'):
        dir_path = dir_path + '/'

    script={}
    #Iterate through list of filenames
    for file in file_names:
        #Get only the name of the script and remove other path information.
        file=file.split('/')[-1]
        file=file.split('.')[0]

        #Import the script file to get the dictionary values
        try:
            script=__import__(file)
        except ValueError:
            print "Error"
        #print script.__maintainer__
        
        #Define output file path
        outf=os.path.join(dir_path,file+'.rst')
        
        
        #This try block attempts to parse the dictionary and if the dictionary
        #is not present, then it will write that information to stdout
        try:
            
            imported_brief_description=script.script_info['brief_description']
            imported_script_description=script.script_info['script_description']

            new_script_description = convert_py_file_to_link(imported_script_description)
            
                
            inputs='\t**[REQUIRED]**\n\t\t\n'

            for i in script.required_options:
                # when no default is provided in the call to make_option,
                # the value of i.default is a tuple -- this try/except
                # handles the diff types that i.default can be
                try:
                    if i.default[0] == 'NO':
                        # i.default is a tuple, so defualt hasn't been
                        # set by the user, and it should therefore be None
                        defaults = None
                    else:
                        # i.default is a string
                        defaults = i.default
                except TypeError:
                    # i.default is not a string or a tuple (e.g., it's an int or None)
                    defaults = i.default
                
                p=re.compile('\%default')
                help_str=p.sub(str(defaults),i.help)
                new_help_str=convert_py_file_to_link(help_str)

                cmd_arg=str(i).replace('--','`-`-').replace('/',' ')
                inputs=inputs+ '\t' + str(cmd_arg) + '\n\t\t' + new_help_str + '\n'
    
            inputs=inputs + '\t\n\t**[OPTIONAL]**\n\t\t\n'
            for i in script.optional_options:
                # when no default is provided in the call to make_option,
                # the value of i.default is a tuple -- this try/except
                # handles the diff types that i.default can be
                try:
                    if i.default[0] == 'NO':
                        
                        # i.default is a tuple, so defualt hasn't been
                        # set by the user, and it should therefore be None
                        defaults = None
                    else:
                        # i.default is a string
                        defaults = i.default
                except TypeError:
                    # i.default is not a string or a tuple (e.g., it's an int or None)
                    defaults = i.default
            
                p=re.compile('\%default')
                help_str=p.sub(str(defaults),i.help)
                new_help_str=convert_py_file_to_link(help_str)
        
                cmd_arg=str(i).replace('--','`-`-').replace('/',', ')
                inputs=inputs+ '\t' + str(cmd_arg) + '\n\t\t' + new_help_str + '\n'

    
            script_examples=''
            for ex in script.script_info['script_usage']:
                if ex[0] <> '':
                    script_examples += '\n**' + ex[0] + '**\n'
                    
                if ex[1] <> '':
                    script_ex=convert_py_file_to_link(ex[1])
                    script_examples += '\n' + script_ex + '\n'
                if ex[2] <>'':
                    script_examples += '\n::\n\n\t' + ex[2] + '\n'
    
            script_out=script.script_info['output_description'].replace('--','`-`-')
            new_script_out=convert_py_file_to_link(script_out)
            

            output_text = rst_text % (file, file, file,         
                                        imported_brief_description,  
                                        new_script_description, file, 
                                        inputs, new_script_out, script_examples)


            ###Write rst file
            f = open(outf, 'w')
            f.write((output_text))
            f.close()
            
            #script.close()
        except AttributeError:
            print "%s: This file does not contain the appropriate dictionary" \
            % (file)


if __name__ == "__main__":
    main()