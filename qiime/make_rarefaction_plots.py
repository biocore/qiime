#!/usr/bin/env python
#file make_rarefaction_plots.py
from __future__ import division
__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2009, QIIME"
__credits__ = ["Meg Pirrung"] 
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
"""

from sys import argv, exit
from random import choice, randrange
from string import strip
from matplotlib.pylab import savefig, clf, gca, gcf, errorbar
import matplotlib.pyplot as plt
import os.path
from os.path import exists, splitext, split
import shutil
from warnings import warn
from itertools import cycle

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

MARKERS = ['*', 'D' , 'H' , 'd' , 'h' , 'o' , 'p' , 's' , 'x']
graphNames = []
sampleIDs = []

def parse_rarefaction_data(lines):
    data = {}
    data['headers'] = []
    data['options'] = []
    data['xaxis'] = []
    data['series'] = {}
    data['error'] = {}
    for l in lines:
        if l.startswith('#'):
            data['headers'].append(l.strip('#').strip())
            continue
        if l.startswith('xaxis'):
            data['xaxis'] = [float(v) for v in l[6:].strip().split('\t')]
            continue
        if l.startswith('>>'):
            data['options'].append(l.strip().strip('>'))
            continue
        if l.startswith('series'):
            data['series'][data['options'][len(data['options'])-1]] = [float(v) for v in l[7:].strip().split('\t')]
            continue
        if l.startswith('error'):
            data['error'][data['options'][len(data['options'])-1]] = [float(v) for v in l[6:].strip().split('\t')]
    # print data
    return data

def save_rarefaction_plots(xaxis, yvals, err, xmax, ymax, ops, mapping_category, \
itype, res, rtype, fpath):
    plt.clf()
    plt.title(rtype + ": " + mapping_category)
    fig  = plt.gcf()
    
    plt.grid(color='gray', linestyle='-')

    ops.sort()
    
    for o in ops:
        l = o
        if len(o) > 20:
            l = l[:20] + '...'
        # plt.errorbar(xaxis[:len(yaxis[o])], yaxis[o], \
        #         yerr=err[o][:len(yaxis[o])], color=colors[o], label=l, \
        #         marker=syms[o], markersize=4)
        
        plt.errorbar(xaxis[:len(yvals[o])], yvals[o], \
        yerr=err[o][:len(yvals[o])], label=l)
        
    c = 1
    if len(ops) > 12:
        c = int(len(ops)/12)

    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, markerscale=.3, ncol=c)
    plt.legend(loc=2, markerscale=.3, ncol=c)
    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.set_xlim((0,xmax))
    # ax.set_ylim((0,ymax))
    ax.set_xlabel('Sequences Per Sample')
    ax.set_ylabel("Rarefaction Measure: " + splitext(split(rtype)[1])[0])
    # plt.show()

    imgpath = fpath +'/'+ mapping_category + '.'+itype
    plt.savefig(imgpath, format=itype, dpi=res)
    imgnmforjs = splitext(split(rtype)[1])[0] + '/' + mapping_category +"."+ itype
    graphNames.append(imgnmforjs)
    
def make_plots(prefs):
    data = {}
    for r in prefs['rarefactions']:
        ymax = prefs['ymax']
        data = parse_rarefaction_data(prefs['rarefactions'][r])
        
        file_path = os.path.join(prefs['output_path'], splitext(split(data['headers'][0])[1])[0])
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        
        if ymax == 0:
            ymax = max([max(v) for v in data['series'].values()])
        xmax = data['xaxis'][len(data['xaxis'])-1]
        
        save_rarefaction_plots(data['xaxis'], data['series'], data['error'], xmax, ymax, data['options'], data['headers'][1], prefs['imagetype'], prefs['resolution'], data['headers'][0], file_path)
    
def make_output_files(prefs, qiime_dir):
    # print 'hurf'
    open(prefs['output_path'] + "/graphNames.txt",'w').writelines([f +\
            '\n' for f in graphNames])
    os.makedirs(prefs['output_path']+"/js")
    os.makedirs(prefs['output_path']+"/css")
    shutil.copyfile(qiime_dir+"/qiime/support_files/html_templates/" + \
    "rarefaction_plots.html", prefs['output_path']+\
    "/rarefaction_plots.html")
    shutil.copyfile(qiime_dir+"/qiime/support_files/js/rarefaction_plots.js", \
    prefs['output_path']+"/js/rarefaction_plots.js")
    shutil.copyfile(qiime_dir+"/qiime/support_files/js/jquery.js", \
    prefs['output_path']+"/js/jquery.js")
    shutil.copyfile(qiime_dir+"/qiime/support_files/js/jquery."+\
    "dataTables.min.js", prefs['output_path']+\
    "/js/jquery.dataTables.min.js")
    shutil.copyfile(qiime_dir+"/qiime/support_files/css/rarefaction_plots.css",\
    prefs['output_path']+"/css/rarefaction_plots.css")
    shutil.copyfile(qiime_dir+"/qiime/support_files/images/qiime_header.png", \
    prefs['output_path']+"/qiime_header.png")