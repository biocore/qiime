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

from string import strip
from matplotlib.pylab import savefig, clf, gca, gcf, errorbar
import matplotlib.pyplot as plt
import os.path
from os.path import exists, splitext, split
import shutil
from itertools import cycle
from qiime.colors import Color, natsort, iter_color_groups


def save_rarefaction_plots(xaxis, yvals, err, xmax, ymax, ops, \
mapping_category, itype, res, rtype, data_colors, colors, fpath,\
 graphNames, background_color):
    plt.clf()
    plt.title((splitext(split(rtype)[1])[0]) + ": " + mapping_category)
    fig  = plt.gcf()
    plt.grid(color='gray', linestyle='-')

    ops = natsort(ops)

    for o in ops:
        l = o
        if len(o) > 20:
            l = l[:20] + '...'
        
        if data_colors != None:
            try:
                plt.errorbar(xaxis[:len(yvals[o])], yvals[o], \
                yerr=err[o][:len(yvals[o])], label=l, color=data_colors[colors[o]].toHex())
            except(KeyError):
                continue
        else:
            plt.errorbar(xaxis[:len(yvals[o])], yvals[o], \
            yerr=err[o][:len(yvals[o])], label=l, color=colors[o])

    c = 1
    if len(ops) > 12:
        c = int(len(ops)/12)

    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, markerscale=.3, ncol=c)
    plt.legend(loc=2, markerscale=.3, ncol=c)
    ax = plt.gca()
    ax.set_axis_bgcolor(background_color)
    ax.set_axisbelow(True)
    ax.set_xlim((0,xmax))
    ax.set_ylim((0,ymax))
    ax.set_xlabel('Sequences Per Sample')
    ax.set_ylabel("Rarefaction Measure: " + splitext(split(rtype)[1])[0])

    imgpath = os.path.join(fpath,mapping_category) + '.'+itype
    plt.savefig(imgpath, format=itype, dpi=res)

    imgnmforjs = splitext(split(rtype)[1])[0] + "/" + mapping_category +"."+\
    itype
    graphNames.append(imgnmforjs)
    return graphNames
    
def make_plots(options):
    data = {}
    graphNames = []
    
    color_prefs,data,background_color,label_color = options['prefs']
    groups_and_colors=iter_color_groups(data['map'],color_prefs)
    groups_and_colors=list(groups_and_colors)
    
    if options['colorby']:
        categs = options['colorby'].split(',')
        options['categories'] = list(set(categs).intersection(set(data['map'][0])))
    
    for r in options['rarefactions']:
        ymax = options['ymax']
        data = options['rarefactions'][r]

        file_path = os.path.join(options['output_path'], \
        splitext(split(data['headers'][0])[1])[0])
        if not os.path.exists(file_path):
            os.makedirs(file_path)

        if ymax == 0:
            ymax = max([max(v) for v in data['series'].values()]) + \
            max([max(v) for v in data['error'].values()])
        xmax = data['xaxis'][len(data['xaxis'])-1]

        if(options['prefs_path'] != None):
            for i in range(len(groups_and_colors)):
                try:
                    groups_and_colors[i].index(r.split('.')[0])
                    break
                except(ValueError):
                    continue
            
            if i == len(groups_and_colors)-1:
                try:
                    groups_and_colors[i].index(r.split('.')[0])
                except(ValueError):
                    continue
                
            groups=groups_and_colors[i][1]
            colors=groups_and_colors[i][2]
            data_colors=groups_and_colors[i][3]

            categories = [k for k in groups]
            
            save_rarefaction_plots(data['xaxis'], data['series'], data['error'], \
            xmax, ymax, categories, data['headers'][1], options['imagetype'], \
            options['resolution'], data['headers'][0], data_colors, colors, file_path, graphNames, background_color)
        
        else:
            if options['colorby']:
                if data['headers'][1] not in options['categories']:
                    continue
                    
            save_rarefaction_plots(data['xaxis'], data['series'], data['error'], \
            xmax, ymax, data['options'], data['headers'][1], options['imagetype'], \
            options['resolution'], data['headers'][0], None, data['color'], file_path, graphNames, background_color)
    return graphNames

def make_output_files(options, qiime_dir, graphNames):
    # print 'hurf'
    os.makedirs(options['output_path']+"/support_files")
    open(options['output_path'] + "/support_files/graphNames.txt",'w').\
    writelines([f +\
    '\n' for f in graphNames])

    os.makedirs(options['output_path']+"/support_files/js")
    os.makedirs(options['output_path']+"/support_files/css")

    shutil.copyfile(qiime_dir+"/qiime/support_files/html_templates/" + \
    "rarefaction_plots.html", options['output_path']+\
    "/rarefaction_plots.html")
    shutil.copyfile(qiime_dir+"/qiime/support_files/js/rarefaction_plots.js", \
    options['output_path']+"/support_files/js/rarefaction_plots.js")
    shutil.copyfile(qiime_dir+"/qiime/support_files/js/jquery.js", \
    options['output_path']+"/support_files/js/jquery.js")
    shutil.copyfile(qiime_dir+"/qiime/support_files/js/jquery."+\
    "dataTables.min.js", options['output_path']+\
    "/support_files/js/jquery.dataTables.min.js")
    shutil.copyfile(qiime_dir+\
    "/qiime/support_files/css/rarefaction_plots.css",\
    options['output_path']+"/support_files/css/rarefaction_plots.css")
    shutil.copyfile(qiime_dir+"/qiime/support_files/images/qiime_header.png", \
    options['output_path']+"/support_files/qiime_header.png")