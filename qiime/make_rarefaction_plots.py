#!/usr/bin/env python
#file make_rarefaction_plots.py
#from __future__ import division
__author__ = "Meg Pirrung and Jesse Stombaugh"
__copyright__ = "Copyright 2009, QIIME"
__credits__ = ["Meg Pirrung","Jesse Stombaugh"] 
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
Matplotlib
Numpy

Example 1: Create rarefaction plot from mapping file and rarefaction data
Usage: python make_rarefaction_plots.py -m mappingfile.txt -r rare1.txt,rare2.txt -p pref1,pref2 -i .png

"""

from sys import argv, exit
from random import choice, randrange
from time import strftime
from qiime import parse
from numpy import array, transpose, random, mean, std
from string import strip
from matplotlib.pylab import savefig, clf, gca
import matplotlib.pyplot as plt
import os.path
from optparse import OptionParser
from os.path import exists, splitext, split
import shutil

COLOUR = ['b', #  : blue
            'r', #  : green
            'g', #  : red
            'y', #  : cyan
            'm', #  : magenta
            'c', #  : yellow
            'k' #  : black 
            ]
err = []

Infinity = 1e10000
NaN = Infinity / Infinity

def is_nan(x):
    return type(x) is float and x != x

def is_finite(x):
    return x != Infinity

def parse_rarefaction(lines):
    col_headers = None
    result = []
    row_headers = []
    for line in lines:
        if line[0] == '#': continue
        if line[0] == '\t': #is header
            col_headers = map(strip, line.split('\t')[1:])
        else:
            entries = line.split('\t')
            try:
                result.append(map(float, entries[1:]))
            except(ValueError):
                result.append(map(strip, map(str, entries[1:])))
                
            row_headers.append(entries[0])

    rare_mat_raw = array(result)
    rare_mat_min = [rare_mat_raw[x][2:] for x in range(0,len(rare_mat_raw))]
    seqs_per_samp = [rare_mat_raw[x][0] for x in range(0,len(rare_mat_raw))]
    sampleIDs = col_headers[2:]
    rare_mat_trans = transpose(array(rare_mat_min)).tolist()
    return rare_mat_trans, seqs_per_samp, sampleIDs
    
def ave_seqs_per_sample(matrix, seqs_per_samp, sampleIDs):
    ave_ser = {}
    for i in range(0,len(sampleIDs)):
        curr = seqs_per_samp[0]
        temp_ser = []
        s = 0
        n = 0
        for j in range(0, len(seqs_per_samp)):
            next = seqs_per_samp[j]
            if curr != next:
                temp_ser.append(s/n)
                s = 0
                n = 0
                curr = next
            else:
                try:
                    s = s + float(matrix[i][j])
                except:
                    s = s + 0
            n = n + 1
        temp_ser.append(s/n)
        ave_ser[sampleIDs[i]] = temp_ser
    #print ave_ser
    return ave_ser

def is_max_category_ops(mapping, mapping_category):
    header = mapping[0][0]
    map_min = mapping[0][1:]
    num_samples = len(map_min)
    index = header.index(mapping_category)
    seen = set()
    
    for m in map_min:
        seen.update([m[index]])
        
    return (len(seen) == num_samples), len(seen)

def make_error_series(rare_mat, sampleIDs, mapping, mapping_category):
    err_ser = dict()
    collapsed_ser = dict()
    header = mapping[0][0]
    map_min = mapping[0][1:]
    mapping_dict = dict()
    for m in map_min:
        mapping_dict[m[0]] = m[1:]
        
    seen = set()
    pre_err = dict()
    op_index = header.index(mapping_category)
    
    for k in mapping_dict.keys():
        try:
            test = rare_mat[k]
        except(KeyError):
            #err.append("SampleID " + k + " found in mapping but not in rarefaction file.\n")
            continue
        op = mapping_dict[k][op_index-1]
        if op not in seen:
            seen.update([op])
            pre_err[op] = []
        pre_err[op].append(rare_mat[k])
    
    #print seen
    ops = [o for o in seen]
    cols = dict() #[COLOUR[i%len(COLOUR)] for i in range(0,len(ops))]
    for i in range(0,len(ops)):
        cols[ops[i]] = COLOUR[i%len(COLOUR)];
        
    for o in ops:
        ao = array(pre_err[o])
        m = mean(ao, 0)
        collapsed_ser[o] = m.tolist()
        s = std(ao, 0)
        err_ser[o] = s.tolist()
    #print ops
    return collapsed_ser, err_ser, ops, cols

def plot_rarefaction_noave(rare_mat, xaxisvals, sampleIDs, mapping, mapping_category):
    xaxis = [float(x) for x in set(xaxisvals)]
    xaxis.sort()
    plt.figure()

    for k in rare_mat.keys():
        rare_mat[k] = [float(v) for v in rare_mat[k] if v != 0]

    for k in rare_mat.keys():
        plt.plot(xaxis[:len(rare_mat[k])], rare_mat[k])

    plt.grid(color='gray', linestyle='-')
    #plt.legend(title=mapping_category)
    ax = plt.gca()
    ax.set_xlabel('Sequences Per Sample')
    return plt

def plot_rarefaction(rare_mat, xaxisvals, sampleIDs, mapping, mapping_category):
    xaxis = [float(x) for x in set(xaxisvals)]
    xaxis.sort()
    plt.figure()
    
    yaxis, err, ops, colors = make_error_series(rare_mat, sampleIDs, mapping, mapping_category)
    plt.axes([.1,.1,.6,.8])
    
    for o in ops:
        try:
            l = o
            if len(o) > 10:
                l = l[:10] + '...'
            plt.errorbar(xaxis, yaxis[o], yerr=err[o], color=colors[o], label=l)
        except(ValueError):
            print mapping_category
            print o
            #print xaxis
            #print yaxis[o]
    plt.grid(color='gray', linestyle='-')
    plt.legend(title=mapping_category, loc=(1.02,.1), markerscale=.5, ncol=int(len(ops)/10)+1)
    ax = plt.gca()
    ax.set_xlabel('Sequences Per Sample')
    return plt

def save_plot(plot, filenm, rtype, title):
    plot.title(title)
    ax = plot.gca()
    ax.set_ylabel(rtype.split('.')[0])
    plot.savefig(filenm +'.png', format='png', dpi=75)

def _make_cmd_parser():
    parser = OptionParser(usage="Usage: make_rarefaction_plots2.py -m <mapping file> \
    -r <rarefaction data file> -p <command line preferences OR preferences file> \
    -i <extension type for image output> -o <output path>")
    
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")
    parser.add_option('-m', '--map', \
        help='name of mapping file [default: %default]')
    parser.add_option('-r', '--rarefaction', \
        help='name of rarefaction file')
    parser.add_option('-p', '--prefs', \
        help='name of columns to make rarefaction graphs of, comma delimited no spaces. Use\
                \'ALL\' command to make graphs of all metadata columns.')
    parser.add_option('-i', '--imagetype', \
        help='extension for image type choose from (.jpg, .gif, .png, .svg, .pdf). DEFAULT (.png)')
    #parser.add_option('-p', '--prefs', \
    #    help='name of preferences file')
    parser.add_option('-o', '--dir_path',\
        help='directory prefix for all analyses [default: %default]',default='')
    options, args = parser.parse_args()
    return options

def get_map(options, data):
    """Opens and returns mapping file for parsing"""
    try:
        #print options.map
        data['map'] = open(options.map, 'U').readlines()
        #print data['map']
        return data['map']
    except (TypeError, IOError):
        print 'Mapping file required for this analysis'
        exit(0)

def get_rarefactions(options, data):
    """Parses and then tries to open rarefaction files to make sure they exist"""
    try:
        #print options.rarefaction
        rarenames = options.rarefaction.split(',')
        rares = dict()
        for r in rarenames:
            rares[r] = open(options.rarefaction, 'U').readlines()
        return rares
    except (TypeError, IOError):
        print 'Rarefaction file required for this analysis'
        exit(0)

def get_prefs(options, data):
    """Opens and returns prefs file for parsing"""
    try:
        if options.prefs.split(',')[0] == 'ALL':
            return options.prefs.split(',')[0]
        ps = options.prefs.split(',')
        for p in ps:
            data['map'][0][0].index(p)
        return ps
    except (ValueError, TypeError, IOError):
        print "Supplied prefs are not found in mapping file, please check spelling and syntax."
        exit(0)

def get_img_extension(options, data):
    """Gets type of extension to save images as."""
    imgtypes = ['.jpg','.gif','.png','.svg','.pdf']
    try:    
        if options.imagetype == None or options.imagetype not in imgtypes:
            print "Extension not supplied or supplied extension not supported, using .png"
            data['imagetype'] = '.png'
        else:
            data['imagetype'] = options.imagetype
        return data['imagetype']
    except (TypeError, IOError):
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
    dir_path = dir_path + 'rarefaction_graphs'

    alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
    alphabet += alphabet.lower()
    alphabet += "01234567890"
    file_path=__file__.split('/')
    qiime_dir = _get_script_dir(argv[0])
    data_file_path=''.join([choice(alphabet) for i in range(10)])
    data_file_path=strftime("%Y_%m_%d_%H_%M_%S")+data_file_path
    data_file_dir_path = dir_path+data_file_path
    data = {}
    
    data['map'] = parse.parse_map(get_map(options,data), return_header=True, strip_quotes=True)
    data['rarefactions'] = get_rarefactions(options,data)
    data['prefs'] = get_prefs(options, data)
    data['output_path'] = data_file_dir_path
    data['imagetype'] = get_img_extension(options, data)
    '''
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
    '''
    return data
    
def make_plots(data):
    for r in data['rarefactions'].keys():
        file_path = data['output_path']+'/'+ r.split('.')[0]
        os.makedirs(file_path)
        try:
            rare_mat_trans, seqs_per_samp, sampleIDs = parse_rarefaction(data['rarefactions'][r])
            rare_mat_ave = ave_seqs_per_sample(rare_mat_trans, seqs_per_samp, sampleIDs)
    
            if data['prefs'] == 'ALL':
                for p in data['map'][0][0]: #headerline
                    is_max, l = is_max_category_ops(data['map'], p)
                    if l == 1:
                        #print "Category \'" + p + "\' only has one option, rarefaction graph was not created."
                        continue
                    if is_max:
                        pr = plot_rarefaction_noave(rare_mat_ave, seqs_per_samp, sampleIDs, data['map'], p)
                    else:
                        pr = plot_rarefaction(rare_mat_ave, seqs_per_samp, sampleIDs, data['map'], p)
                    filenm = file_path + '/'+ p
                    save_plot(pr, filenm, r, r.split('.')[0] +'_'+ p)
                    clf()
            else:
                for p in data['prefs']:
                    is_max, l = is_max_category_ops(data['map'], p)
                    if l == 1:
                        #print "Category \'" + p + "\' only has one option, rarefaction graph was not created."
                        continue
                    if is_max:
                        pr = plot_rarefaction_noave(rare_mat_ave, seqs_per_samp, sampleIDs, data['map'], p)
                    else:
                        pr = plot_rarefaction(rare_mat_ave, seqs_per_samp, sampleIDs, data['map'], p)
                    filenm = file_path + '/'+ p
                    save_plot(pr, filenm, r, r.split('.')[0] +'_'+ p)
                    clf()
        except():
            os.removedirs(file_path)
            print "Error:", sys.exc_info()[0]
    
if __name__ == '__main__':
    from sys import argv, exit
    options = _make_cmd_parser()
    file_data = _process_prefs(options)
    make_plots(file_data)
