#!/usr/bin/env python
# File created on 14 Jul 2014
from __future__ import division

from os.path import exists, splitext, join, isdir
from os import makedirs, listdir, remove
from qiime.util import RExecutor
from biom import load_table

__author__ = "Sophie Weiss"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Sophie Weiss", "Joseph Paulson"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Sophie Weiss"
__email__ = "sophie.sjw@gmail.com"



def normalize_CSS(input_path, out_path, output_CSS_statistics):
    """performs metagenomeSeq's CSS normalization on a single raw abundance OTU matrix
    """
    base_fname, ext = splitext(input_path)
    json_infile = base_fname+'_json.biom'
    open(str(json_infile),'w').write(load_table(input_path).to_json('forR'))

    if output_CSS_statistics:
        base_fname, ext = splitext(out_path)
        output_CSS_statistics = base_fname+'_CSS_statistics.txt'
    
    run_CSS(json_infile, out_path, output_CSS_statistics=output_CSS_statistics)
    remove(json_infile)


def multiple_file_normalize_CSS(input_dir, output_dir, output_CSS_statistics):
    """performs metagenomeSeq's CSS normalization on a directory of raw abundance OTU matrices
    """
    if not exists(output_dir):
        makedirs(output_dir)
    file_names = [fname for fname in listdir(input_dir) if not (fname.startswith('.')\
        or isdir(fname))]

    for fname in file_names:
        base_fname, ext = splitext(fname)
        print base_fname
        original_fname = base_fname+'.biom'
        json_fname = base_fname+'_json.biom'
        hdf5_infile = join(input_dir, original_fname)
        json_infile = join(input_dir, json_fname)
        open(str(json_infile),'w').write(load_table(hdf5_infile).to_json('forR'))
        outfile = join(output_dir, 'CSS_'+base_fname+'.biom')
        if output_CSS_statistics:
            output_CSS_statistics = join(output_dir, 'CSS_statistics_'+base_fname+'.txt')
        run_CSS(json_infile, outfile, output_CSS_statistics=output_CSS_statistics)
        remove(json_infile)

def run_CSS(input_path, out_path, output_CSS_statistics, HALT_EXEC=False):
    """Run metagenomeSeq's fitZIG algorithm through Rscript
    """
    # set options
    if output_CSS_statistics==False:
        command_args = ['-i %s -o %s' % (input_path, out_path)]
    else:
        command_args = ['-i %s -o %s -s %s' % (input_path, out_path, output_CSS_statistics)]
    # instantiate the object
    rsl = RExecutor()
    # run the app
    app_result = rsl(command_args=command_args, script_name='CSS.r')

    return app_result


def normalize_DESeq(input_path, out_path, DESeq_negatives_to_zero):
    """performs DESeqVS normalization on a single raw abundance OTU matrix
    """
    base_fname, ext = splitext(out_path)
    json_infile = base_fname+'_json.biom'
    open(str(json_infile),'w').write(load_table(input_path).to_json('forR'))

    run_DESeq(json_infile, out_path, DESeq_negatives_to_zero)
    remove(json_infile)

def multiple_file_normalize_DESeq(input_dir, output_dir, DESeq_negatives_to_zero):
    """performs DESeqVS normalization on a directory of raw abundance OTU matrices
    """
    if not exists(output_dir):
        makedirs(output_dir)
    file_names = [fname for fname in listdir(input_dir) if not (fname.startswith('.')\
        or isdir(fname))]

    for fname in file_names:
        base_fname, ext = splitext(fname)
        original_fname = base_fname+'.biom'
        json_fname = base_fname+'_json.biom'
        hdf5_infile = join(input_dir, original_fname)
        json_infile = join(input_dir, json_fname)
        open(str(json_infile),'w').write(load_table(hdf5_infile).to_json('forR'))
        outfile = join(output_dir, 'DESeqVS_'+base_fname+'.biom')

        run_DESeq(json_infile, outfile, DESeq_negatives_to_zero)
        remove(json_infile)

def run_DESeq(input_path, out_path, DESeq_negatives_to_zero, HALT_EXEC=False):
    """Run metagenomeSeq's fitZIG algorithm through Rscript
    """
    # set options
    if DESeq_negatives_to_zero:
        command_args = ['-i %s -o %s -z %s' % (input_path, out_path, DESeq_negatives_to_zero)]
    else:
        command_args = ['-i %s -o %s' % (input_path, out_path)]
    # instantiate the object
    rsl = RExecutor()
    # run the app
    app_result = rsl(command_args=command_args, script_name='DESeq2.r')

    return app_result

def algorithm_list():
    """ returns list of normalization algorithms from qiime.normalize_table
    """
    return ['CSS', 'DESeq']

if __name__ == "__main__":
    main()
