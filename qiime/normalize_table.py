#!/usr/bin/env python
# File created on 14 Jul 2014
from __future__ import division

from os.path import exists, splitext, join, isdir
from os import makedirs, listdir, remove

import tempfile
from biom import load_table
from qiime.util import RExecutor
from qiime.util import get_qiime_temp_dir


__author__ = "Sophie Weiss"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Sophie Weiss", "Joseph Paulson"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Sophie Weiss"
__email__ = "sophie.sjw@gmail.com"



def normalize_CSS(input_path, out_path, output_CSS_statistics):
    """performs metagenomeSeq's CSS normalization on a single raw abundance OTU matrix
    """
    tmp_bt = load_table(input_path) 

    if output_CSS_statistics:
        base_fname, ext = splitext(out_path)
        output_CSS_statistics = base_fname+'_CSS_statistics.txt'

    with tempfile.NamedTemporaryFile(dir=get_qiime_temp_dir(),
                                     prefix='QIIME-normalize-table-temp-table-',
                                     suffix='.biom') as temp_fh:
        temp_fh.write(tmp_bt.to_json('forR'))
        temp_fh.flush()
        run_CSS(temp_fh.name, out_path, output_CSS_statistics=output_CSS_statistics)


def multiple_file_normalize_CSS(input_dir, output_dir, output_CSS_statistics):
    """performs metagenomeSeq's CSS normalization on a directory of raw abundance OTU matrices
    """
    if not exists(output_dir):
        makedirs(output_dir)
    file_names = [fname for fname in listdir(input_dir) if not (fname.startswith('.')\
        or isdir(fname))]

    for fname in file_names:
        base_fname, ext = splitext(fname)
        original_fname = base_fname+'.biom'
        hdf5_infile = join(input_dir, original_fname)
        tmp_bt = load_table(hdf5_infile) 
        outfile = join(output_dir, 'CSS_'+base_fname+'.biom')
        if output_CSS_statistics:
            output_CSS_statistics = join(output_dir, 'CSS_statistics_'+base_fname+'.txt')

        with tempfile.NamedTemporaryFile(dir=get_qiime_temp_dir(),
                                         prefix='QIIME-normalize-table-temp-table-',
                                         suffix='.biom') as temp_fh:
            temp_fh.write(tmp_bt.to_json('forR'))
            temp_fh.flush()
            run_CSS(temp_fh.name, outfile, output_CSS_statistics=output_CSS_statistics)

def run_CSS(input_path, out_path, output_CSS_statistics):
    """Run metagenomeSeq's CSS algorithm through Rscript
    """
    # set options
    if not output_CSS_statistics:
        command_args = ['-i %s -o %s' % (input_path, out_path)]
    else:
        command_args = ['-i %s -o %s -s %s' % (input_path, out_path, output_CSS_statistics)]
    # instantiate the object
    rsl = RExecutor(TmpDir=get_qiime_temp_dir())
    # run the app
    app_result = rsl(command_args=command_args, script_name='CSS.r')

    return app_result


def normalize_DESeq2(input_path, out_path, DESeq_negatives_to_zero):
    """performs DESeq2VS normalization on a single raw abundance OTU matrix
    """
    tmp_bt = load_table(input_path) 
    with tempfile.NamedTemporaryFile(dir=get_qiime_temp_dir(),
                                     prefix='QIIME-normalize-table-temp-table-',
                                     suffix='.biom') as temp_fh:
        temp_fh.write(tmp_bt.to_json('forR'))
        temp_fh.flush()
        run_DESeq2(temp_fh.name, out_path, DESeq_negatives_to_zero)

def multiple_file_normalize_DESeq2(input_dir, output_dir, DESeq_negatives_to_zero):
    """performs DESeq2VS normalization on a directory of raw abundance OTU matrices
    """
    if not exists(output_dir):
        makedirs(output_dir)
    file_names = [fname for fname in listdir(input_dir) if not (fname.startswith('.')\
        or isdir(fname))]

    for fname in file_names:
        base_fname, ext = splitext(fname)
        original_fname = base_fname+'.biom'
        hdf5_infile = join(input_dir, original_fname)
        tmp_bt = load_table(hdf5_infile) 
        outfile = join(output_dir, 'DESeq2_'+base_fname+'.biom')

        with tempfile.NamedTemporaryFile(dir=get_qiime_temp_dir(),
                                         prefix='QIIME-normalize-table-temp-table-',
                                         suffix='.biom') as temp_fh:
            temp_fh.write(tmp_bt.to_json('forR'))
            temp_fh.flush()
            run_DESeq2(temp_fh.name, outfile, DESeq_negatives_to_zero)

def run_DESeq2(input_path, out_path, DESeq_negatives_to_zero):
    """Run DESeq2's variance stabilization algorithm through Rscript
    """
    # set options
    if DESeq_negatives_to_zero:
        command_args = ['-i %s -o %s -z %s' % (input_path, out_path, DESeq_negatives_to_zero)]
    else:
        command_args = ['-i %s -o %s' % (input_path, out_path)]
    # instantiate the object
    rsl = RExecutor(TmpDir=get_qiime_temp_dir())
    # run the app
    app_result = rsl(command_args=command_args, script_name='DESeq2.r')

    return app_result

def algorithm_list():
    """ returns list of normalization algorithms from qiime.normalize_table
    """
    return ['CSS', 'DESeq2']

if __name__ == "__main__":
    main()
