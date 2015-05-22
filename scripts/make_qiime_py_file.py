#!/usr/bin/env python
# make_qiime_py_file.py

"""
This is a script which will add headers and footers to new python files
and make them executable.
"""

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from sys import exit
from os import popen
from os.path import exists
from time import strftime
from optparse import OptionParser

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Create python file"""
script_info['script_description'] = ("This is a script which will add headers "
"and footers to new python files and make them executable.")
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example usage:""",
     """Create a new script:""",
     '%prog -s -a "Greg Caporaso" -e gregcaporaso@gmail.com -o my_script.py'))
script_info['script_usage'].append(
    ('',
     'Create a new test file:',
     '%prog -t -a "Greg Caporaso" -e gregcaporaso@gmail.com -o my_test.py'))
script_info['script_usage'].append(
    ('',
     'Create a basic file (e.g., for library code):',
     '%prog -a "Greg Caporaso" -e gregcaporaso@gmail.com -o my_lib.py'))
script_info['script_usage_output_to_remove'] = [
    'my_lib.py',
    'my_script.py',
    'my_test.py']
script_info['output_description'] = ("The results of this script is either a "
                                     "python script, test, or library file, "
                                     "depending on the input parameters.")
script_info['required_options'] = [
    options_lookup['output_fp']
]

script_info['optional_options'] = [
    make_option('-s', '--script', action='store_true', default=False,
                help=("Pass if creating a script to include option parsing "
                      "framework [default:%default].")),
    make_option('-t', '--test', action='store_true', default=False,
                help=("Pass if creating a unit test file to include relevant "
                      "information [default:%default].")),
    make_option('-a', '--author_name', default='AUTHOR_NAME',
                help=("The script author's (probably you) name to be included " 
                      "the header variables. This will typically need to be "
                      "enclosed in quotes to handle spaces. "
                      "[default:%default]")),
    make_option('-e', '--author_email', default='AUTHOR_EMAIL',
                help=("The script author's (probably you) e-mail address to "
                      "be included the header variables. [default:%default]")),
    make_option('-c', '--copyright',
                default='Copyright 2014, The QIIME Project',
                help=("The copyright information to be included in "
                      "the header variables. [default:%default]"))
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    header_block = """#!/usr/bin/env python
# File created on %s
from __future__ import division

__author__ = "AUTHOR_NAME"
__copyright__ = "COPYRIGHT"
__credits__ = ["AUTHOR_NAME"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "AUTHOR_NAME"
__email__ = "AUTHOR_EMAIL"
""" % strftime('%d %b %Y')

    script_block = """
from qiime.util import parse_command_line_parameters, make_option

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
# Members of the tuple in script_usage are (title, description, example call)
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
    # Example required option
    #make_option('-i', '--input_fp', type='existing_filepath',
    #            help='the input filepath')
]
script_info['optional_options'] = [
    # Example optional option
    #make_option('-o', '--output_dir', type="new_dirpath",
    #            help='the output directory [default: %default]')
]
script_info['version'] = __version__"""

    test_block = """
from shutil import rmtree
from os.path import exists, join
from tempfile import mkstemp
from os import close

from cogent.util.unit_test import TestCase, main
from skbio.util import remove_files, create_dir

from qiime.util import get_qiime_temp_dir
from qiime.test import initiate_timeout, disable_timeout

class NAMETests(TestCase):

    def setUp(self):
        \"\"\" \"\"\"
        self.files_to_remove = []
        self.dirs_to_remove = []

        # Create example output directory
        tmp_dir = get_qiime_temp_dir()
        fd, self.test_out = mkstemp(dir=tmp_dir,
                                   prefix='qiime_parallel_tests_',
                                   suffix='')
        close(fd)
        self.dirs_to_remove.append(self.test_out)
        create_dir(self.test_out)

        # Create example input file
        fd, self.inseqs1_fp = mkstemp(dir=self.test_out,
                                     prefix='qiime_inseqs',
                                     suffix='.fasta')
        close(fd)
        inseqs1_f = open(self.inseqs1_fp,'w')
        inseqs1_f.write(inseqs1)
        inseqs1_f.close()
        self.files_to_remove.append(self.inseqs1_fp)

        # Define number of seconds a test can run for before timing out
        # and failing
        initiate_timeout(60)

    def tearDown(self):
        \"\"\" \"\"\"
        disable_timeout()
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

inseqs1 = \"\"\">example input here
ACGT
\"\"\""""

    if opts.test and opts.script:
        option_parser.error('-s and -t cannot both be passed: file must be a '
                            'test or a script, or neither one.')

    script = opts.script
    test = opts.test
    output_fp = opts.output_fp

    # Check to see if the file which was requested to be created
    # already exists -- if it does, print a message and exit
    if exists(output_fp):
        option_parser.error("The file name you requested already exists. "
                            "Delete existing file and rerun script if it "
                            "should be overwritten. Otherwise change the file "
                            "name (-o). Creating no files and exiting...")

    # Create the header data
    header_block = header_block.replace('AUTHOR_NAME', opts.author_name)
    header_block = header_block.replace('AUTHOR_EMAIL', opts.author_email)
    header_block = header_block.replace('COPYRIGHT', opts.copyright)
    lines = [header_block]

    if test:
        lines.append(test_block)
        lines += ['', '', 'if __name__ == "__main__":', '    main()']
    elif script:
        lines.append(script_block)
        lines += ['', '', 'def main():',
                  '    option_parser, opts, args = '
                  'parse_command_line_parameters(**script_info)',
                  '',
                  'if __name__ == "__main__":',
                  '    main()']

    # Open the new file for writing and write it.
    f = open(output_fp, 'w')
    f.write('\n'.join(lines))
    f.close()

    if test or script:
        # Change the permissions on the new file to make it executable
        chmod_string = ' '.join(['chmod 755', output_fp])
        popen(chmod_string)

if __name__ == "__main__":
    main()
