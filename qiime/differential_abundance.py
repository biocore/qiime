#!/usr/bin/env python
# File created on 14 Jul 2014
from __future__ import division

from os.path import exists, splitext, join, isdir
from os import makedirs, listdir, remove, path

from qiime.parse import parse_mapping_file_to_dict
from qiime.util import RExecutor
from biom import load_table


__author__ = "Sophie Weiss"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Sophie Weiss", "Joseph Paulson"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Sophie Weiss"
__email__ = "sophie.sjw@gmail.com"



def DA_fitZIG(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2):
   """perform metagenomeSeq's Zero Inflated Gaussian (ZIG) OTU differential abundance testing"""
   base_fname, ext = splitext(input_path)
   json_infile = base_fname+'_json.biom'
   tmp_bt = load_table(input_path)
   tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
   tmp_bt.add_metadata(tmp_pmf, 'sample')
   open(str(json_infile),'w').write(tmp_bt.to_json('forR'))

   run_fitZIG(json_infile, out_path, mapping_category, subcategory_1, subcategory_2)
   remove(json_infile)


def multiple_file_DA_fitZIG(input_dir, output_dir, mapping_fp, mapping_category, subcategory_1, subcategory_2):
    """perform metagenomeSeq's Zero Inflated Gaussian (ZIG) OTU differential abundance test on a directory of raw abundance OTU matrices
    """
    if not exists(output_dir):
        makedirs(output_dir)
    file_names = listdir(input_dir)
    file_names = [fname for fname in file_names if not (fname.startswith('.')\
        or isdir(fname))]

    for fname in file_names:
        base_fname, ext = splitext(fname)
        original_fname = base_fname+'.biom'
        json_fname = base_fname+'_json.biom'
        hdf5_infile = join(input_dir, original_fname)
        json_infile = join(input_dir, json_fname)
        tmp_bt = load_table(hdf5_infile)
        #add metadata for R
        tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
        tmp_bt.add_metadata(tmp_pmf, 'sample')
        #make temporary json biom version - R currently does not have hdf5
        open(str(json_infile),'w').write(tmp_bt.to_json('forR'))
        outfile = join(output_dir, 'fitZIG_DA_'+base_fname+'.txt')

        run_fitZIG(json_infile, outfile, mapping_category, subcategory_1, subcategory_2)
        #delete temporary json formatted OTU table  
        remove(json_infile)


def run_fitZIG(input_path, out_path, mapping_category, subcategory_1, subcategory_2, HALT_EXEC=False):
    """Run metagenomeSeq's fitZIG algorithm through Rscript
    """
    # set options
    command_args = ['-i %s -o %s -c %s -x %s -y %s' % (input_path, out_path, mapping_category, subcategory_1, subcategory_2)]
    # instantiate the object
    rsl = RExecutor()
    # run the app
    app_result = rsl(command_args=command_args, script_name='fitZIG.r')

    return app_result


def DA_DESeq2(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots):
    """perform DESeq2 negative binomial differential abundance test on a raw abundance OTU matrix
    """
    base_fname, ext = splitext(input_path)
    json_infile = base_fname+'_json.biom'
    tmp_bt = load_table(input_path)
    #add metadata for R
    tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
    tmp_bt.add_metadata(tmp_pmf, 'sample')
    #make temporary json biom version - R currently does not have hdf5
    open(str(json_infile),'w').write(tmp_bt.to_json('forR'))
    base_fname, ext = splitext(out_path)
    outfile_diagnostic = join(base_fname+'_diagnostic_plots.pdf') 

    run_DESeq2(json_infile, out_path, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots, outfile_diagnostic)
    remove(json_infile)


def multiple_file_DA_DESeq2(input_dir, output_dir, mapping_fp, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots):
    """perform DESeq2 negative binomial differential abundance test on a directory of raw abundance OTU matrices
    """
    if not exists(output_dir):
        makedirs(output_dir)
    file_names = listdir(input_dir)
    file_names = [fname for fname in file_names if not (fname.startswith('.')\
        or isdir(fname))]

    for fname in file_names:
        base_fname, ext = splitext(fname)
        print base_fname
        original_fname = base_fname+'.biom'
        json_fname = base_fname+'_json.biom'
        hdf5_infile = join(input_dir, original_fname)
        json_infile = join(input_dir, json_fname)
        tmp_bt = load_table(hdf5_infile)
        tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
        tmp_bt.add_metadata(tmp_pmf, 'sample')
        open(str(json_infile),'w').write(tmp_bt.to_json('forR'))
        outfile = join(output_dir, 'DESeq2_DA_'+base_fname+'.txt') 
        outfile_diagnostic = join(output_dir, 'DESeq2_diagnostic_plots_'+base_fname+'.pdf') 

        run_DESeq2(json_infile, outfile, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots, outfile_diagnostic)
        remove(json_infile)


def run_DESeq2(input_path, out_path, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots, outfile_diagnostic, HALT_EXEC=False):
    """Run metagenomeSeq's fitZIG algorithm through Rscript
    """
    # set options
    if DESeq2_diagnostic_plots==True:
        command_args = ['-i %s -o %s -c %s -x %s -y %s -d %s -e %s' % (input_path, out_path, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots, outfile_diagnostic)]
    else:
        command_args = ['-i %s -o %s -c %s -x %s -y %s' % (input_path, out_path, mapping_category, subcategory_1, subcategory_2)]
    # instantiate the object
    rsl = RExecutor()
    # run the app
    app_result = rsl(command_args=command_args, script_name='DESeq2_nbinom.r')

    return app_result


def algorithm_list():
    """ returns list of differential abundance detection algorithms from qiime.differential_abundance
    """
    return ['metagenomeSeq_fitZIG', 'DESeq2_nbinom']