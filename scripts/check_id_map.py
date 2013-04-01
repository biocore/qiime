#!/usr/bin/env python
# File created on 27 Oct 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from qiime.command import cl_main
from qiime.cl_interfaces import CheckIdMap
from sys import argv

cmd = CheckIdMap()
script_info = cmd.getScriptInfo()
if __name__ == "__main__":
    cl_main(cmd,argv)