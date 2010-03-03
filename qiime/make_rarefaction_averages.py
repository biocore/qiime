#!/usr/bin/env python
#file make_rarefaction_averages.py
from __future__ import division
__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2009, QIIME"
__credits__ = ["Meg Pirrung"] 
__license__ = "GPL"
__version__ = "0.91"
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

from sys import argv, exit
from random import choice, randrange
from time import strftime
from qiime import parse, util
from numpy import array, transpose, random, mean, std, arange
from string import strip
import os.path
from os.path import exists, splitext, split
import shutil
from warnings import warn
from itertools import cycle

def parse_rarefaction(lines):
    """Function for parsing rarefaction files specifically for use in
    make_rarefaction_averages.py"""
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
                temp = []
                for x in entries[1:]:
                    if x.strip() != 'n/a':
                        temp.append(float(x.strip()))
                    else:
                        temp.append(0.0)
                result.append(temp)
                
            row_headers.append(entries[0])
    rare_mat_raw = array(result)
    rare_mat_min = [rare_mat_raw[x][2:] for x in range(0,len(rare_mat_raw))]
    seqs_per_samp = [rare_mat_raw[x][0] for x in range(0,len(rare_mat_raw))]
    sampleIDs = col_headers[2:]
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
    header = mapping[0][0]
    map_min = mapping[0][1:]
    num_samples = len(map_min)
    index = header.index(mapping_category)
    seen = set()
    
    for m in map_min:
        seen.update([m[index]])
        
    return (len(seen) == num_samples), len(seen)

def make_error_series(rtype, rare_mat, sampleIDs, mapping, mapping_category):
    """Create mean and error bar series for the supplied mapping category"""
    err_ser = dict()
    collapsed_ser = dict()
    header = mapping[0][0]
    map_min = mapping[0][1:]
    mapping_dict = dict()
    notfound = []
    for m in map_min:
        mapping_dict[m[0]] = m[1:]
        
    seen = set()
    pre_err = {}
    category_index = header.index(mapping_category)
    
    mapkeys = set(mapping_dict.keys())
    rarekeys = set(rare_mat.keys())
    
    intsc = mapkeys.intersection(rarekeys)
    diff = mapkeys.difference(rarekeys)
    
    if len(intsc) == 0:
        print "Error: None of the sampleIDs found in the mapping file are"+\
        " in the rarefaction file %s, check to make sure these files"%rtype+\
        " correspond."
        exit(0)
    elif len(diff) > 0:
        warnstr = "SampleIDs %s found in mapping but not in "%list(diff)+\
        "rarefaction file %s."%rtype
        warn(warnstr)
    
    for k in intsc:
        test = rare_mat[k]
        op = mapping_dict[k][category_index-1]
        if op not in seen:
            seen.update([op])
            pre_err[op] = []
        pre_err[op].append(rare_mat[k])

    
    ops = list(seen)
    
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
mapping_category, rare_type, fpath):
    yaxis, err, ops = make_error_series(rare_type, rare_mat, \
    sampleIDs, mapping, mapping_category)
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
    for o in ops:
        lines.append(">> " + o + '\n')
        lines.append('series ')
        line = ''
        for v in yaxis[o]:
            line += str(v) + '\t'
        line += '\n'
        lines.append(line)
        lines.append('error ')
        line = ''
        for e in err[o]:
            line += str(e) + '\t'
        line += '\n'
        lines.append(line)
    open(fpath+'/'+mapping_category+'.txt','w').writelines(lines)

def make_averages(prefs):
    rarelines = []

    for r in prefs['rarefactions']:
        file_path = os.path.join(prefs['output_path'], \
        splitext(split(r)[1])[0])
        os.makedirs(file_path)
        
        rare_mat_trans = prefs['rarefactions'][r][0]
        seqs_per_samp = prefs['rarefactions'][r][1]
        sampleIDs = prefs['rarefactions'][r][2]

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
            
        for p in prefs['categories']:
            save_rarefaction_data(rare_mat_ave, xaxisvals, xmax, \
            sampleIDs, prefs['map'], p, r, file_path)
        
    tablelines = ['#SampleIDs\n']
    tablelines.extend([s + '\n' for s in sampleIDs])
    tablelines.extend(rarelines)
    open(prefs['output_path'] + "/rarefactionTable.txt",'w').writelines(tablelines)