#!/usr/bin/env python
# File created on 22 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from qiime.util import make_option
from qiime.util import parse_command_line_parameters, get_options_lookup

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Deprecation warning: split_libraries_illumina.py is deprecated. You should use split_libraries_fastq.py. See the 'Processing Illumina Data tutorial at: http://qiime.org/tutorials/processing_illumina_data.html."
script_info['script_description'] = "Deprecation warning: split_libraries_illumina.py is deprecated. You should use split_libraries_fastq.py. See the 'Processing Illumina Data tutorial at: http://qiime.org/tutorials/processing_illumina_data.html."
script_info['script_usage'] = []
script_info['output_description']= ""
script_info['required_options'] = []
script_info['optional_options'] = []
script_info['version'] = __version__


def main():
    print "Deprecation warning: split_libraries_illumina.py is deprecated. You should use split_libraries_fastq.py. See the 'Processing Illumina Data tutorial at: http://qiime.org/tutorials/processing_illumina_data.html."
    
if __name__ == "__main__":
    main()
    