#!/usr/bin/env python
#file make_rarefaction_averages.py
from __future__ import division
__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Meg Pirrung"] 
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Meg Pirrung"
__email__ = "meg.pirrung@colorado.edu"
__status__ = "Release"

"""
Author: Meg Pirrung (meg.pirrung@colorado.edu) 
Status: Prototype

Requirements:
Python 2.5
Matplotlib
Numpy
"""

from sys import exit
from random import choice, randrange
from qiime.parse import parse_rarefaction
from qiime.colors import iter_color_groups
from numpy import array, transpose, random, mean, std, arange
from string import strip
import os.path
from os.path import exists, splitext, split
from warnings import warn

def get_rarefaction_data(rarefaction_data, col_headers):
    rare_mat_raw = array(rarefaction_data)
    rare_mat_min = [rare_mat_raw[x][2:] for x in range(0,len(rare_mat_raw))]
    seqs_per_samp = [rare_mat_raw[x][0] for x in range(0,len(rare_mat_raw))]
    sampleIDs = col_headers[3:]
    rare_mat_trans = transpose(array(rare_mat_min)).tolist()
    return rare_mat_trans, seqs_per_samp, sampleIDs

def ave_seqs_per_sample(matrix, seqs_per_samp, sampleIDs):
    """Calculate the average for each sampleID across each number of \
    seqs/sample"""
    ave_ser = {}
    temp_dict = {}
    for i,sid in enumerate(sampleIDs):
        temp_dict[sid] = {}
        for j,seq in enumerate(seqs_per_samp):
            try:
                temp_dict[sid][seq].append(matrix[i][j])
            except(KeyError):
                temp_dict[sid][seq] = []
                temp_dict[sid][seq].append(matrix[i][j])

    for sid in sampleIDs:
        ave_ser[sid] = []
        keys = temp_dict[sid].keys()
        keys.sort()
        for k in keys:
            ave_ser[sid].append(mean(array(temp_dict[sid][k]),0))
    return ave_ser

def is_max_category_ops(mapping, mapping_category):
    """Count how many unique values there are for the supplied mapping \
    category and return true if all values are unique"""
    header = mapping[1]
    map_min = mapping[0]
    num_samples = len(map_min)
    index = header.index(mapping_category)
    seen = set()
    
    for m in map_min:
        seen.update([m[index]])
        
    return (len(seen) == num_samples), len(seen)

def make_error_series(rtype, rare_mat, groups, mapping_category):
    """Create mean and error bar series for the supplied mapping category"""
    err_ser = dict()
    collapsed_ser = dict()
        
    seen = set()
    pre_err = {}
    
    ops = [k for k in groups]
    
    
    # mapkeys = set([])
    #     rarekeys = set(rare_mat.keys())
    #     
    #     intsc = mapkeys.intersection(rarekeys)
    #     diff = mapkeys.difference(rarekeys)
    #     
    #     if len(intsc) == 0:
    #         print "Error: None of the sampleIDs found in the mapping file are"+\
    #         " in the rarefaction file %s, check to make sure these files"%rtype+\
    #         " correspond."
    #         exit(0)
    #     elif len(diff) > 0:
    #         warnstr = "SampleIDs %s found in mapping but not in "%list(diff)+\
    #         "rarefaction file %s."%rtype
    #         warn(warnstr)
    
    # for k in intsc:
    #     test = rare_mat[k]
    #     op = mapping_dict[k][header.index(mapping_category)]
    #     if op not in seen:
    #         seen.update([op])
    #         pre_err[op] = []
    #     pre_err[op].append(rare_mat[k])

    notfound = []

    for o in ops:
        pre_err[o] = []
        for samID in groups[o]:
            try:
                pre_err[o].append(rare_mat[samID])
            except(KeyError):
                notfound.append(samID)
    
    # print "The following sampleIDs were found in the mapping file but not in \
    # the rarefaction file: "
    # print notfound
    
    for o in ops:
        min_len = 1000 #1e10000
        for series in pre_err[o]:
            series = [float(v) for v in series if v != 0]
            if len(series) < min_len:
                min_len = len(series)
        
        pre_err[o] = [x[:min_len] for x in pre_err[o]]
    
    for o in ops:
        opsarray = array(pre_err[o])
        mn = mean(opsarray, 0)
        collapsed_ser[o] = mn.tolist()
        stddev = std(opsarray, 0)
        err_ser[o] = stddev.tolist()
        
    return collapsed_ser, err_ser, ops

def get_overall_averages(rare_mat, sampleIDs):
    """Make series of averages of all values of seqs/sample for each \
    sampleID"""
    overall_ave = dict();
    for s in sampleIDs:
        overall_ave[s] = mean(array(rare_mat[s]))
    return overall_ave

def save_rarefaction_data(rare_mat, xaxis, xmax, sampleIDs, mapping, \
mapping_category, colors, rare_type, fpath, data_colors, groups):
    yaxis, err, ops = make_error_series(rare_type, rare_mat, \
    groups, mapping_category)
    # print colors.keys()
    lines = []
    lines.append("# "+rare_type+'\n')
    lines.append("# "+mapping_category+'\n')
    line = ''
    line += 'xaxis: '
    for v in xaxis:
        line += str(v) + '\t'
    line += '\n'
    lines.append(line)
    lines.append('xmax: '+str(xmax)+'\n')
    
    for o in colors.keys():
        lines.append(">> " + o + '\n')
        if colors != None:
            try:
                lines.append("color " + data_colors[colors[o]].toHex() + '\n')
            except(KeyError):
                print 'wtf'

        lines.append('series ')
        line = ''
        try:
            for v in yaxis[o]:
                line += str(v) + '\t'
        except(TypeError):
            line += str(yaxis[o])
        line += '\n'
        lines.append(line)
        
        lines.append('error ')
        line = ''
        try:
            for e in err[o]:
                line += str(e) + '\t'
        except(TypeError):
            line += str(err[o])
        line += '\n'
        lines.append(line)
    open(fpath+'/'+mapping_category+'.txt','w').writelines(lines)

def make_averages(options):
    rarelines = []
    
    color_prefs,data,background_color,label_color = options['prefs']

    groups_and_colors=iter_color_groups(data['map'],color_prefs)
    groups_and_colors=list(groups_and_colors)

    if options['colorby']:
        categs = options['colorby'].split(',')
        options['categories'] = list(set(categs).intersection(set(data['map'][0])))

    for r in options['rarefactions']:
        file_path = os.path.join(options['output_path'], \
        splitext(split(r)[1])[0])
        os.makedirs(file_path)
        
        col_headers, comments, rarefaction_fn, rarefaction_data = \
        options['rarefactions'][r]
        rare_mat_trans, seqs_per_samp, sampleIDs = \
        get_rarefaction_data(rarefaction_data, col_headers)

        xaxisvals = [float(x) for x in set(seqs_per_samp)]
        xaxisvals.sort()

        rare_mat_ave = ave_seqs_per_sample(rare_mat_trans, seqs_per_samp, \
        sampleIDs)
        xmax = max(xaxisvals) + (xaxisvals[len(xaxisvals)-1] - \
        xaxisvals[len(xaxisvals)-2])
        overall_average = get_overall_averages(rare_mat_ave, sampleIDs)

        rarelines.append("#" + r + '\n')
        
        for s in sampleIDs:
            rarelines.append('%f'%overall_average[s] + '\n')
        
        if options['prefs']:
            for i in range(len(groups_and_colors)):
                labelname=groups_and_colors[i][0]
                groups=groups_and_colors[i][1]
                colors=groups_and_colors[i][2]
                data_colors=groups_and_colors[i][3]
                data_color_order=groups_and_colors[i][4]
            
                save_rarefaction_data(rare_mat_ave, xaxisvals, xmax, \
                sampleIDs, data['map'], labelname, colors, r, \
                file_path, data_colors, groups)
        else:
            for p in options['cagtegories']:
                save_rarefaction_data(rare_mat_ave, xaxisvals, xmax, \
                sampleIDs, data['map'], p, None, r, file_path, None)
        
    tablelines = ['#SampleIDs\n']
    tablelines.extend([s + '\n' for s in sampleIDs])
    tablelines.extend(rarelines)
    open(options['output_path'] + "/rarefactionTable.txt",'w').\
    writelines(tablelines)