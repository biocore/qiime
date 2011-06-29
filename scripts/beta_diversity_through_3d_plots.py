#!/usr/bin/env python
# File created on 04 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"


deprecation_message = """Usage: beta_diversity_through_3d_plots.py

beta_diversity_through_3d_plots.py has been renamed in beta_diversity_through_plots.py, incorporating several new features.

beta_diversity_through_3d_plots.py was removed from QIIME in version 1.2.1-dev.

No analysis has been performed and no output has been written."""

script_info={}
script_info['brief_description']=""""""
script_info['script_description']= deprecation_message
script_info['script_usage']=[]
script_info['output_description']=""
script_info['required_options']=[]
script_info['optional_options']=[]

script_info['version'] = __version__


def main():
    print deprecation_message

if __name__ == "__main__":
    main()