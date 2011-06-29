#!/usr/bin/env python
# File created on 10 Feb 2011.
from __future__ import division

__author__ = "Jesse Stombauhg"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"


deprecation_message = """Usage: make_pie_charts.py

make_pie_charts.py has been incorporated in plot_taxa_summary.py, where the user
can generate pie, area and bar taxa summaries.

make_pie_charts.py was removed from QIIME in version 1.2.0-dev.

No analysis performed and no output written."""

script_info={}
script_info['brief_description']=""""""
script_info['script_description']="""Usage: make_pie_charts.py

make_pie_charts.py has been incorporated in plot_taxa_summary.py, where the user
can generate pie, area and bar taxa summaries.

make_pie_charts.py was removed from QIIME in version 1.2.0-dev.

No analysis performed and no output written."""
script_info['script_usage']=[]
script_info['output_description']=""
script_info['required_options']=[]
script_info['optional_options']=[]

script_info['version'] = __version__

def main():
    print deprecation_message

if __name__ == "__main__":
    main()