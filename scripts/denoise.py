#!/usr/bin/env python
# File created on 10 Feb 2011.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"


deprecation_message = """Usage: denoise.py

denoise.py has been replaced by denoise_wrapper.py in QIIME 1.2.1-dev. Call denoise_wrapper.py -h for usage information.

No analysis performed and no output written."""

script_info={}
script_info['brief_description']=""""""
script_info['script_description']="""Usage: denoise.py

denoise.py has been replaced by denoise_wrapper.py in QIIME 1.2.1-dev. Call denoise_wrapper.py -h for usage information.


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