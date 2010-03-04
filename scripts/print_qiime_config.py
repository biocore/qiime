#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder","Dan Knights"]
__license__ = "GPL"
__version__ = "0.92"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Release"

from qiime.util import load_qiime_config, parse_command_line_parameters
from optparse import make_option

script_info = {}
script_info['brief_description']= """Print out the qiime config settings."""
script_info['script_description'] = """A simple scripts that prints out the qiime config settings."""
script_info['script_usage']=[]
script_info['script_usage'].append(("Example 1","""Print qiime config settings:""","""print_qiime_config.py"""))
script_info['output_description'] = """This prints the qiime_config to stdout."""
script_info['version'] = __version__
script_info['help_on_no_arguments'] = False

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    qiime_config = load_qiime_config()
    for key,value in  qiime_config.items():
        print "%20s:\t%s"%(key,value)

if __name__ == "__main__":
    main()
