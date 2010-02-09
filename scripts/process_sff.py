#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.process_sff import prep_sffs_in_dir

script_description = """Converts directory of sff files into fasta and qual files.

Requires that 454's off-instrument apps (sffinfo, sfffile) are on your path.
"""

script_usage = """Convert all the sffs in directory sffs to fasta and qual:

process_sff -i sffs
"""

required_options = [\
 make_option('-i','--input_dir',help='the input directory of sffs')
]

optional_options = [\
]




def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    prep_sffs_in_dir(opts.input_dir)

if __name__ == "__main__":
    main()
