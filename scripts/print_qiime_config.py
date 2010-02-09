#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Pre-release"

from qiime.util import load_qiime_config, parse_command_line_parameters
from optparse import make_option

script_description = \
"""A simple scripts that prints out the qiime config settings."""

script_usage = """To print qiime config settings:
print_qiime_config """

required_options = []

optional_options = []

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options,
      help_on_no_arguments=False)

    qiime_config = load_qiime_config()
    for key,value in  qiime_config.items():
        print "%20s:\t%s"%(key,value)

if __name__ == "__main__":
    main()
