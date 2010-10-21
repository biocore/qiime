#!/usr/bin/env python
#file process_sff.py
from cogent.util.misc import app_path
from cogent.app.util import ApplicationNotFoundError
from os import listdir, system
from os.path import splitext, join, isfile, isdir,split
"""Converts directory of sff files into fasta and qual files.

Requires that 454's off-instrument apps (sffinfo, sfffile) are on your path.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight","Greg Caporaso", "Jesse Stombaugh"] 
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

def check_sfffile():
    """Raise error if sfffile is not in $PATH """
    if not app_path('sfffile'):
        raise ApplicationNotFoundError,\
        "sfffile is not in $PATH. Is it installed? Have you added it to $PATH?"

def convert_Ti_to_FLX(filename,output_pathname):
    """Converts Titanium SFF to FLX length reads."""
    check_sfffile()
    system('sfffile -flx -o %s %s' % (output_pathname,filename))

def make_flow_txt(filename,output_pathname):
    """Makes flowgram file from sff file."""
    check_sffinfo()
    system('sffinfo %s > %s.txt' % (filename, output_pathname))

def make_fna(filename,output_pathname):
    """Makes fna file from sff file."""
    check_sffinfo()
    system('sffinfo -s %s > %s.fna' % (filename, output_pathname))

def make_qual(filename,output_pathname):
    """Makes qual file from sff file."""
    check_sffinfo()
    system('sffinfo -q %s > %s.qual' % (filename, output_pathname))

def prep_sffs_in_dir(pathname,make_flowgram, output_pathname,convert_to_flx):
    """Converts all sffs in dir to fasta/qual."""
    check_sffinfo()
    
    # this statement is for passing a single sff file
    if isfile(pathname):
        name=split(pathname)[-1]
        pathname=split(pathname)[:-1][0]
        # if this an sff do the following
        if name.endswith('.sff'):
            new_pathname=join(pathname,name)
            
            # if the user wants to convert to FLX len we need to write a new
            # sff file
            if convert_to_flx:
                check_sfffile()
                FLX_fname=join(output_pathname,splitext(name)[0]+'_FLX.sff')
                convert_Ti_to_FLX(join(pathname,name),FLX_fname)
                new_pathname=FLX_fname
                name=splitext(name)[0]+'_FLX.sff'
            
            # write fna and qual files
            make_fna(new_pathname,join(output_pathname,splitext(name)[0]))
            make_qual(new_pathname,join(output_pathname,splitext(name)[0]))
            
            #write flow file
            if make_flowgram:
                make_flow_txt(new_pathname,join(output_pathname,splitext(name)[0]))
                
    # this statement is for passing a directory containing sff files
    elif isdir(pathname):
        for name in listdir(pathname):
            
            if name.endswith('.sff'):
                new_pathname=join(pathname,name)
                # if the user wants to convert to FLX len we need to write a new
                # sff file
                if convert_to_flx:
                    check_sfffile()
                    FLX_fname=join(output_pathname,splitext(name)[0]+'_FLX.sff')
                    convert_Ti_to_FLX(join(pathname,name),FLX_fname)
                    new_pathname=FLX_fname
                    name=splitext(name)[0]+'_FLX.sff'
                
                # write fna and qual files
                make_fna(new_pathname,join(output_pathname,splitext(name)[0]))
                make_qual(new_pathname,join(output_pathname,splitext(name)[0]))
                #write flow file
                if make_flowgram:
                    make_flow_txt(new_pathname,join(output_pathname,splitext(name)[0]))
    else:
        raise OSError, "The path '%s' is not valid!" % pathname
        
