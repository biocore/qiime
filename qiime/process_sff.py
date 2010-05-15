#!/usr/bin/env python
#file process_sff.py
from cogent.util.misc import app_path
from cogent.app.util import ApplicationNotFoundError
from os import listdir, system
from os.path import splitext, join
"""Converts directory of sff files into fasta and qual files.

Requires that 454's off-instrument apps (sffinfo, sfffile) are on your path.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight","Greg Caporaso"] 
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

def check_sffinfo():
    """Raise error if sffinfo is not in $PATH """
    if not app_path('sffinfo'):
        raise ApplicationNotFoundError,\
         "sffinfo is not in $PATH. Is it installed? Have you added it to $PATH?"

def make_fna(filename):
    """Makes fna file from sff file."""
    check_sffinfo()
    system('sffinfo -s %s > %s.fna' % (filename, splitext(filename)[0]))

def make_qual(filename):
    """Makes qual file from sff file."""
    check_sffinfo()
    system('sffinfo -q %s > %s.qual' % (filename, splitext(filename)[0]))

def prep_sffs_in_dir(pathname):
    """Converts all sffs in dir to fasta/qual."""
    check_sffinfo()
    for name in listdir(pathname):
        if name.endswith('.sff'):
            make_fna(join(pathname,name))
            make_qual(join(pathname,name))
