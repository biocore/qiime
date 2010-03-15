#!/usr/bin/env python
# File created on 15 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

from os.path import join, abspath
from random import choice
from datetime import datetime
from os import makedirs
 
def get_random_directory_name(suppress_mkdir=False,\
    timestamp_pattern='%Y%m%d%H%M%S',\
    rand_length=20,\
    output_dir=None,\
    prefix='',
    suffix=''):
    """Build a random directory name and create the directory 
    
        suppress_mkdir: only build the directory name, don't
         create the directory (default: False)
        timestamp_pattern: string passed to strftime() to generate
         the timestamp (pass '' to suppress the timestamp)
        rand_length: length of random string of characters
        output_dir: the directory which should contain the 
         random directory
        prefix: prefix for directory name
        suffix: suffix for directory name
    
    """
    output_dir = output_dir or './'
    # Define a set of characters to be used in the random directory name
    chars = "abcdefghigklmnopqrstuvwxyz"
    picks = chars + chars.upper() + "0123456790"
    
    # Get a time stamp
    timestamp = datetime.now().strftime(timestamp_pattern)
        
    # Construct the directory name
    dirname = '%s%s%s%s' % (prefix,timestamp,\
                        ''.join([choice(picks) for i in range(rand_length)]),\
                        suffix)
    dirpath = abspath(join(output_dir,dirname))
    
    # Make the directory
    if not suppress_mkdir:
        try:
            makedirs(dirpath)
        except OSError:
            raise OSError,\
             "Cannot make directory %s. Do you have write access?" % dirpath
             
    # Return the path to the directory
    return dirpath
