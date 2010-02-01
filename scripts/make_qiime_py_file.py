#!/usr/bin/env python
# make_qiime_py_file.py

"""
This is a script which will add headers and footers to new python files
and make them executable. 
"""

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

from sys import exit
from os import popen
from os.path import exists
from time import strftime
from optparse import OptionParser

usage_str = """usage: %prog [options] {-o OUTPUT_FP}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
Create a new script:
 python %prog -s -a "Greg Caporaso" -e gregcaporaso@gmail.com -o my_script.py
 
Create a new test file:
 python %prog -t -a "Greg Caporaso" -e gregcaporaso@gmail.com -o my_test.py
 
Create a basic file (e.g., for library code):
 python %prog -a "Greg Caporaso" -e gregcaporaso@gmail.com -o my_lib.py
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)
     
    parser.add_option('-o','--output_fp',
     help="The file path to create [REQUIRED]")

    parser.add_option('-s','--script',action='store_true',\
     help="Pass if creating a script to include option parsing"+\
     " framework [default:%default].")

    parser.add_option('-t','--test',action='store_true',\
     help="Pass if creating a unit test file to include relevant"+\
     " information [default:%default].")
     
    parser.add_option('-a','--author_name',
     help="The script author's (probably you) name to be included"+\
     " the header variables. This will typically need to be enclosed "+\
     " in quotes to handle spaces. [default:%default]")
     
    parser.add_option('-e','--author_email',
     help="The script author's (probably you) e-mail address to be included"+\
     " the header variables. [default:%default]")
     
    parser.add_option('-c','--copyright',
     help="The copyright information to be included in"+\
     " the header variables. [default:%default]")

    parser.set_defaults(script=False,test=False,
     author_name='AUTHOR_NAME',\
     author_email='AUTHOR_EMAIL',\
     copyright='Copyright 2010, The QIIME project')
    
    opts, args = parser.parse_args()
    required_options = ['output_fp']
    
    if opts.test and opts.script:
        parser.error('-s and -t cannot both be passed: file must be a test'+\
                     ' or a script, or neither one.')
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 
     
    return opts, args

script_block = """
from optparse import OptionParser

usage_str = \"\"\"usage: %prog [options] {required options}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
\"\"\"

def parse_command_line_parameters():
    \"\"\" Parses command line arguments \"\"\"
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    # A binary 'verbose' flag
    parser.add_option('-v','--verbose',action='store_true',\\
        dest='verbose',help='Print information during execution -- '+\\
        'useful for debugging [default: %default]')

    # An example REQUIRED option
    #parser.add_option('-i','--input_dir',\\
    #     help='the input directory [REQUIRED]')
    
    # An example option
    #parser.add_option('-o','--output_dir',\\
    #     help='the output directory [default: %default]')

    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False)

    opts,args = parser.parse_args()
    # list of required options (e.g., required_options = ['input_dir'])
    required_options = []
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 
            
    return opts, args

"""

header_block =\
"""#!/usr/bin/env python
# File created on %s
from __future__ import division

__author__ = "AUTHOR_NAME"
__copyright__ = "COPYRIGHT"
__credits__ = ["AUTHOR_NAME"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "AUTHOR_NAME"
__email__ = "AUTHOR_EMAIL"
__status__ = "Pre-release"
 
""" % strftime('%d %b %Y')

if __name__ == "__main__":
    
    opts,args = parse_command_line_parameters()
    script = opts.script
    test = opts.test
    output_fp = opts.output_fp

    # Check to see if the file which was requested to be created
    # already exists -- if it does, print a message and exit
    if exists(output_fp):
        print '\n'.join(["The file name you requested already exists.",\
            " Delete extant file and rerun script if it should be overwritten.",\
            " Otherwise change the file name (-o).",\
            "Creating no files and exiting..."])
        exit(1) 


    # Create the header data
    header_block = header_block.replace('AUTHOR_NAME',opts.author_name)
    header_block = header_block.replace('AUTHOR_EMAIL',opts.author_email)
    header_block = header_block.replace('COPYRIGHT',opts.copyright)
    lines = [header_block]

    # If this is a test file, add the requiste import statement
    if test:
        lines.append('from cogent.util.unit_test import TestCase, main')
    elif script:
        lines.append(script_block)

    # Create the footer data
    # File will be executable
    lines += ['','','','if __name__ == "__main__":']
    if test:
        # Run unittest.main() if test file
        lines.append('    main()')
    elif script:
        lines.append('    opts,args = parse_command_line_parameters()')
        lines.append('    verbose = opts.verbose')
    else:
        # Running the file does nothing by default if not a test file
        lines.append('    pass')

    # Open the new file for writing and write it.
    f = open(output_fp,'w')
    f.write('\n'.join(lines))
    f.close()

    # Change the permissions on the new file to make it executable
    chmod_string = ' '.join(['chmod 755',output_fp])
    popen(chmod_string)
