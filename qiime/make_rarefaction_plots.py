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
Usage: python make_rarefaction_plots.py -m mappingfile.txt -r rare1.txt,rare2.txt
or
python make_rarefaction_plots.py -m mappingfile.txt -r rare1.txt,rare2.txt -p SampleID -i png -d 150

"""

from sys import argv, exit
from random import choice, randrange
from time import strftime
from qiime import parse
from numpy import array, transpose, random, mean, std, arange
from string import strip
from matplotlib.pylab import savefig, clf, gca, gcf, errorbar
import matplotlib.pyplot as plt
import os.path
from optparse import OptionParser
from os.path import exists, splitext, split
import shutil

COLOUR_GRAD = ['#9933cc', #purple
        '#3333cc', #blue
        '#6699cc', #bluetint
        '#666666', #gray
        '#9966cc', #lilac
        '#009999', #cyan
        '#66cc99', #lightgreen
        '#3399cc', #skyblue
        '#00cc66', #sea green
        '#33cc33', #green
        '#99cc00', #yellowgreen
        '#cccc00', #yellow
        '#cc6600', #orange
        '#cc6633', #peach
        '#cc3300', #dark orange
        '#cc6666', #coral
        '#cc0066', #hot pink
        '#cc0000', #red
        ]

COLOUR = ['#9933cc', #purple
            '#99cc00', #yellowgreen
            '#3399cc', #skyblue
            '#66cc99', #lightgreen
            '#ffff00', #yellow
            '#cc0066', #hot pink
            '#cc6600', #orange
            '#00cc66', #sea green
            '#6699cc', #bluetint
            '#cc6633', #peach
            '#33cc33', #green
            '#9966cc', #lilac
            '#3333cc', #blue
            '#cccc00', #gold
            '#cc0000', #red
            '#cc6666', #coral
            '#666666', #gray
            '#009999', #cyan
]

#MARKERS = ['+' , '*' , ',' , '.' , '1' , '2' , '3' , '4' , '<' , '>' , 'D' , 'H' , '^' , '_' , 'd' , 'h' , 'o' , 'p' , 's' , 'v' , 'x' , '|']
MARKERS = ['*', 'D' , 'H' , 'd' , 'h' , 'o' , 'p' , 's' , 'x']

'''COLOUR = ['b', #  : blue
            'r', #  : green
            'g', #  : red
            'y', #  : cyan
            'm', #  : magenta
            'c', #  : yellow
            'k' #  : black 
            ]'''
err = []
graphNames = []
sampleIDs = []
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
                # fix so that float works on all vals except n/a and leave n/as in
                # then transpose matrix, get rid of na's, average
                result.append(map(str, entries[1:]))
                
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
        s = 0 #sum
        n = 0 #iterator
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
    cols = dict()
    syms = dict()
    for i in range(0,len(ops)):
        cols[ops[i]] = COLOUR[i%len(COLOUR)]
        syms[ops[i]] = MARKERS[i%len(MARKERS)]
    
    for o in ops:
        min_len = 100
        for s in pre_err[o]:
            s = [float(v) for v in s if v != 'n/a' and v != 0]
            if len(s) < min_len:
                min_len = len(s)
        
        pre_err[o] = [x[:min_len] for x in pre_err[o]]
            
    for o in ops:
        try:
            ao = array(pre_err[o])
            m = mean(ao, 0)
            collapsed_ser[o] = m.tolist()
            s = std(ao, 0)
            err_ser[o] = s.tolist()
        except(ValueError):
            continue
            print o
        
    return collapsed_ser, err_ser, ops, cols, syms

def get_overall_averages(rare_mat, sampleIDs):
    overall_ave = dict();
    for s in sampleIDs:
        overall_ave[s] = mean(array(rare_mat[s]))
    return overall_ave

def plot_rarefaction_noave(rare_mat, xaxis, sampleIDs, mapping, mapping_category):
    plt.gcf().set_size_inches(8,6)    
    plt.grid(color='gray', linestyle='-')
    
    yseries = []
    for k in rare_mat.keys():
        yseries.append([float(v) for v in rare_mat[k] if v != 'NA' and v != 0])
        
    for s in yseries:
        plt.plot(xaxis[:len(s)], s)

    plt.grid(color='gray', linestyle='-')
    ax = plt.gca()
    ax.set_xlabel('Sequences Per Sample')
    return plt

def plot_rarefaction(rare_mat, xaxis, sampleIDs, mapping, mapping_category):
    plt.gcf().set_size_inches(10,6)   
    plt.grid(color='gray', linestyle='-')
    
    yaxis, err, ops, colors, syms = make_error_series(rare_mat, sampleIDs, mapping, mapping_category)
    ops.sort()
    
    for o in ops:
        try:
            yaxis[o] = [float(v) for v in yaxis[o] if v != 'NA' and v != 0]
            l = o
            if len(o) > 20:
                l = l[:20] + '...'
            plt.errorbar(xaxis[:len(yaxis[o])], yaxis[o], yerr=err[o][:len(yaxis[o])], color=colors[o], label=l, marker=syms[o], markersize=4)
            test = errorbar(xaxis[:len(yaxis[o])], yaxis[o], yerr=err[o][:len(yaxis[o])], color=colors[o], label=l, marker=syms[o], markersize=4)
        except(ValueError):
            print mapping_category
            print o
            #print xaxis
            #print yaxis[o]
    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.set_xlabel('Sequences Per Sample')
    return plt, ops, colors, syms, test

def save_plot(plot, filenm, rtype, title, itype, res, xmax, ymax, ops, cols, syms, line):
    plot.title(title)
    ax = plot.gca()
    ax.set_xlim((0,xmax))
    ax.set_ylim((0,ymax))
    ax.set_ylabel("Rarefaction Type: " + splitext(split(rtype)[1])[0])
    plot.savefig(filenm +'.'+itype, format=itype, dpi=res)
    plt.clf()
    if ops != None:
        c = 1
        if len(ops) > 12:
            c = int(len(ops)/12)
        i = 0
        for o in ops:
            plt.plot([0,0], [0,0], cols[o], label=o, marker=syms[o])
            i += 1
        plt.legend(markerscale=.3, ncol=c)
    plot.savefig(filenm +'_legend.'+itype, format=itype, dpi=res)

def _make_cmd_parser():
    parser = OptionParser(usage="Usage: make_rarefaction_plots.py -m <mapping file> \
    -r <rarefaction data file> -p <command line preferences OR preferences file> \
    -i <extension type for image output> -d <resolution for output image> -o <output path>")

    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")
    parser.add_option('-m', '--map', \
        help='name of mapping file [REQUIRED]')
    parser.add_option('-r', '--rarefaction', \
        help='name of rarefaction file [REQUIRED]')
    parser.add_option('-p', '--prefs', \
        help='name of columns to make rarefaction graphs of, comma delimited no spaces. Use \'ALL\' command to make graphs of all metadata columns. [default: %default]', default='ALL')
    parser.add_option('-i', '--imagetype', \
        help='extension for image type choose from (jpg, gif, png, svg, pdf). [default: %default]', default='png')
    #parser.add_option('-p', '--prefs', \
    #    help='name of preferences file')
    parser.add_option('-d', '--resolution', \
        help='output image resolution, [default: %default]', default="75")
    parser.add_option('-o', '--dir_path',\
        help='directory prefix for all analyses [default: %default]',default='.')
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
            rares[r] = open(r, 'U').readlines()
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
    imgtypes = ['jpg','gif','png','svg','pdf']
    try:    
        if options.imagetype not in imgtypes:
            print "Supplied extension not supported, using .png instead."
            data['imagetype'] = 'png'
        else:
            data['imagetype'] = options.imagetype
        return data['imagetype']
    except (TypeError, IOError):
        return None
        
def get_resolution(options, data):
    """Gets image resolution."""
    try:    
        try:
            data['resolution'] = int(options.resolution)
        except(ValueError):
            print "Inavlid resolution, proceeding with 75dpi."
            data['resolution'] = 75
        return data['resolution']
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
    dir_path = options.dir_path
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
    data['map'][0][0] = [h.strip('#').strip(' ') for h in data['map'][0][0]]
    data['rarefactions'] = get_rarefactions(options,data)
    data['prefs'] = get_prefs(options, data)
    data['output_path'] = data_file_dir_path
    data['imagetype'] = get_img_extension(options, data)
    data['resolution'] = get_resolution(options, data)
    return data, qiime_dir
    
def make_plots(data):
    rarelines = []
    
    for r in data['rarefactions'].keys():
        file_path = data['output_path']+'/'+ splitext(split(r)[1])[0]
        os.makedirs(file_path)
        try:
            rare_mat_trans, seqs_per_samp, sampleIDs = parse_rarefaction(data['rarefactions'][r])
            
            xaxisvals = [float(x) for x in set(seqs_per_samp)]
            xaxisvals.sort()
            
            rare_mat_ave = ave_seqs_per_sample(rare_mat_trans, seqs_per_samp, sampleIDs)
            xmax = max(xaxisvals) + (xaxisvals[len(xaxisvals)-1] - xaxisvals[len(xaxisvals)-2])
            ymax = max([max(s) for s in rare_mat_ave.values()]) + 5
            overall_average = get_overall_averages(rare_mat_ave, sampleIDs)
            rarelines.append("#" + r + '\n')
            for s in sampleIDs:
                rarelines.append('%f'%overall_average[s] + '\n')
    
            if data['prefs'] == 'ALL':
                for p in data['map'][0][0]: #headerline
                    is_max, l = is_max_category_ops(data['map'], p)
                    if l == 1:
                        #print "Category \'" + p + "\' only has one option, rarefaction graph was not created."
                        continue
                    if is_max:
                        pr = plot_rarefaction_noave(rare_mat_ave, xaxisvals, sampleIDs, data['map'], p)
                        ops = None
                        cols = None
                        syms = None
                        line = None
                    else:
                        pr,ops,cols,syms,line = plot_rarefaction(rare_mat_ave, xaxisvals, sampleIDs, data['map'], p)
                    filenm = file_path + '/'+ p
                    graphNames.append(splitext(split(r)[1])[0] + '/'+p+"."+data['imagetype'])
                    save_plot(pr, filenm, r, splitext(split(r)[1])[0] +':'+ p, data['imagetype'], data['resolution'], xmax, ymax, ops, cols, syms, line)
                    plt.clf()
            else:
                for p in data['prefs']:
                    is_max, l = is_max_category_ops(data['map'], p)
                    if l == 1:
                        # print "Category \'" + p + "\' only has one option, rarefaction graph was not created."
                        continue
                    if is_max:
                        pr = plot_rarefaction_noave(rare_mat_ave, xaxisvals, sampleIDs, data['map'], p)
                        ops = None
                        cols = None
                        line = None
                    else:
                        pr,ops,cols,syms,line = plot_rarefaction(rare_mat_ave, xaxisvals, sampleIDs, data['map'], p)
                    filenm = file_path + '/'+ p
                    graphNames.append(splitext(split(r)[1])[0] + '/'+p+"."+data['imagetype'])
                    save_plot(pr, filenm, r, splitext(split(r)[1])[0] +': '+ p, data['imagetype'], data['resolution'], xmax, ymax, ops, cols, syms, line)
                    plt.clf()
        except():
            os.removedirs(file_path)
            print "Error:", sys.exc_info()[0]
            
    tablelines = ['#SampleIDs\n']
    tablelines.extend([s + '\n' for s in sampleIDs])
    tablelines.extend(rarelines)
    return tablelines
    
def make_output_files(data, lines, qiime_dir):
    open(data['output_path'] + "/graphNames.txt",'w').writelines([f +'\n' for f in graphNames])
    open(data['output_path'] + "/rarefactionTable.txt",'w').writelines(lines)
    
    os.makedirs(data['output_path']+"/js")
    os.makedirs(data['output_path']+"/css")
    # open(data['output_path'] + "/rarefaction_plots.html",'w').writelines(open(qiime_dir + "rarefaction_plots.html", "U").readlines())
    # open(data['output_path'] + "/js/rarefaction_plots.js",'w').writelines(open(qiime_dir + "/js/rarefaction_plots.js", "U").readlines())
    # open(data['output_path'] + "/js/jquery.js",'w').writelines(open(qiime_dir + "/js/jquery.js", "U").readlines())
    # open(data['output_path'] + "/js/jquery.dataTables.min.js",'w').writelines(open(qiime_dir + "/js/jquery.dataTables.min.js", "U").readlines())
    # open(data['output_path'] + "/css/rarefaction_plots.css",'w').writelines(open(qiime_dir + "/css/rarefaction_plots.css", "U").readlines())
    shutil.copyfile(qiime_dir+"/rarefaction_plots.html", data['output_path']+"/rarefaction_plots.html")
    shutil.copyfile(qiime_dir+"/js/rarefaction_plots.js", data['output_path']+"/js/rarefaction_plots.js")
    shutil.copyfile(qiime_dir+"/js/jquery.js", data['output_path']+"/js/jquery.js")
    shutil.copyfile(qiime_dir+"/js/jquery.dataTables.min.js", data['output_path']+"/js/jquery.dataTables.min.js")
    shutil.copyfile(qiime_dir+"/css/rarefaction_plots.css", data['output_path']+"/css/rarefaction_plots.css")
    shutil.copyfile(qiime_dir+"/qiime_header.png", data['output_path']+"/qiime_header.png")
    
if __name__ == '__main__':
    from sys import argv, exit
    options = _make_cmd_parser()
    file_data, parent_directory = _process_prefs(options)
    outputlines = make_plots(file_data)
    make_output_files(file_data, outputlines, parent_directory)
