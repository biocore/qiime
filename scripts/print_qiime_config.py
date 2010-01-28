#!/usr/bin/env python 

"""A simple scripts that prints out the qiime config settings."""

from optparse import OptionParser
from qiime.util import load_qiime_config

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

def parse_command_line_parameters(commandline_args=None):
    """returns command-line options"""
    from sys import argv

    usage = """print_qiime_config.py

Prints out the values of all qimme variables stored in the qiime_config file."""

    version = '%prog ' + str(__version__)
    parser = OptionParser(usage=usage, version=version)
    opts,args = parser.parse_args(commandline_args)
    
    return opts,args

def main(commandline_args=None):

    opts, args = parse_command_line_parameters(commandline_args)
    qiime_config = load_qiime_config()
    for key,value in  qiime_config.items():
        print "%20s:\t%s"%(key,value)

if __name__ == "__main__":

    main()
