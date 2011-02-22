#!/usr/bin/env python
# File created on 19 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"


deprecation_message = """Usage: jackknifed_upgma.py

jackknifed_upgma.py has been renamed jackknifed_beta_diversity.py.

jackknifed_upgma.py will be removed from QIIME in version 1.2.x-dev.

No analysis performed and no output written."""

script_info={}
script_info['brief_description']=""""""
script_info['script_description']="""Usage: jackknifed_upgma.py

jackknifed_upgma.py has been renamed jackknifed_beta_diversity.py.

jackknifed_upgma.py will be removed from QIIME in version 1.2.x-dev.

No analysis performed and no output written."""
script_info['script_usage']=[]
script_info['output_description']=""
script_info['required_options']=[]
script_info['optional_options']=[]

script_info['version'] = __version__


print deprecation_message
