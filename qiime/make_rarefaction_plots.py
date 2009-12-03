#!/usr/bin/env python
#file make_rarefaction_plots.py

__author__ = "Meg Pirrung and Jesse Stombaugh"
__copyright__ = "Copyright 2009, the 454 Project" #consider project name
__credits__ = ["Meg Pirrung","Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Meg Pirrung"
__email__ = "meg.pirrung@colorado.edu"
__status__ = "Prototype"

"""
Author: Meg Pirrung (meg.pirrung@colorado.edu) 
Status: Prototype

Requirements:
Python 2.5

Example 1: Create rarefaction plot from mapping file and rarefaction data
Usage: python make_rarefaction_plots.py -m mappingfile.txt -r rare1.txt,rare2.txt


"""

from sys import argv    
from random import choice, randrange
from time import strftime
import os.path
from optparse import OptionParser
from os.path import exists, splitext, split
import shutil

def _make_cmd_parser():
    parser = OptionParser(usage="Usage: rfFileMaker.py -m <mapping file> \
    -r <rarefaction data files, comma separated, no spaces> -o <output path>")
    
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")
    parser.add_option('-m', '--map_fname', \
        help='name of mapping file [default: %default]')
    parser.add_option('-r', '--rarefactions', \
        help='name of rarefaction files comma delimited, no spaces')
    parser.add_option('-o', '--dir_path',\
        help='directory prefix for all analyses [default: %default]',default='')
    options, args = parser.parse_args()
    return options

def get_map(options, data):
    """Tries to open mapping file to make sure it exists"""
    try:
        open(options.map_fname, 'U')
        data['map'] = options.map_fname
        return data['map']
    except (TypeError, IOError):
        print 'Mapping file required for this analysis'
        return None

def get_rarefactions(options, data):
    """Parses and then tries to open rarefaction files to make sure they exist"""
    try:
        #print options.rarefactions
        filenms = options.rarefactions.split(',')
        for f in filenms:
            open(f)
        data['rarefactions'] = filenms
        return data['rarefactions']
    except (TypeError, IOError):
        print 'At least one rarefaction file required for this analysis'
        return None

def _get_script_dir(script_path):
    """Returns directory current script is running in.
    """
    if '/' in script_path:
        script_dir = script_path.rsplit('/',1)[0]+'/'
    else:
        script_dir = './'
    return script_dir

def _process_prefs(options):    
    dir_path = "."
    if dir_path and not dir_path.endswith('/'):
        dir_path = dir_path + '/'
    dir_path = dir_path + 'rarefaction_ui'

    alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
    alphabet += alphabet.lower()
    alphabet += "01234567890"
    
    file_path=__file__.split('/')
    
    qiime_dir = _get_script_dir(argv[0])

    data_file_path=''.join([choice(alphabet) for i in range(10)])
    data_file_path=strftime("%Y_%m_%d_%H_%M_%S")+data_file_path
    data_file_dir_path = dir_path+data_file_path

    #data_file_dir_path=create_dir(data_file_dir_path,'')
    os.mkdir(data_file_dir_path)
    os.mkdir(data_file_dir_path + '/js')
    os.mkdir(data_file_dir_path + '/css')
    
    data = {}
    
    data['map'] = get_map(options,data)
    data['rarefactions'] = get_rarefactions(options,data)
    
    filenms = []
    filenms.append(data['map'])
    for r in data['rarefactions']:
        filenms.append(r.split('/')[len(r.split('/'))-1])
    #print filenms
    
    open(data_file_dir_path + "/dataFilesforJS.txt",'w').writelines([f.split('/')[len(f.split('/'))-1]+'\n' for f in filenms])
    for r in data['rarefactions']:
        fl = open(data_file_dir_path + "/" + r.split('/')[len(r.split('/'))-1],'w')
        fl.writelines([l.replace('\r','\n').strip(' ') for l in open(r).readlines()])
    open(data_file_dir_path + "/rarefaction_plots.html",'w').writelines(open(qiime_dir + "/rarefaction_plots.html").readlines())
    open(data_file_dir_path + "/js/rarefaction_plots.js",'w').writelines(open(qiime_dir + "/js/rarefaction_plots.js").readlines())
    open(data_file_dir_path + "/js/jquery.js",'w').writelines(open(qiime_dir + "/js/jquery.js").readlines())
    open(data_file_dir_path + "/js/jquery.dataTables.min.js",'w').writelines(open(qiime_dir + "/js/jquery.dataTables.min.js").readlines())
    open(data_file_dir_path + "/css/rarefaction_plots.css",'w').writelines(open(qiime_dir + "/css/rarefaction_plots.css").readlines())
    
    shutil.copyfile(qiime_dir+"/qiime_header.png", data_file_dir_path+"/qiime_header.png")

if __name__ == '__main__':
    from sys import argv, exit
    options = _make_cmd_parser()
    _process_prefs(options)