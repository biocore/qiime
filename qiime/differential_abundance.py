#!/usr/bin/env python
# File created on 14 Jul 2014
from __future__ import division

from os.path import exists, splitext, join, isdir
from os import makedirs, listdir, remove, path

import tempfile
from qiime.util import get_qiime_temp_dir

from qiime.parse import parse_mapping_file_to_dict
from qiime.util import RExecutor, MetadataMap
from biom import load_table


__author__ = "Sophie Weiss"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Sophie Weiss", "Joseph Paulson"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Sophie Weiss"
__email__ = "sophie.sjw@gmail.com"


def check_mapping_file_category(loaded_biom, mapping_fp, mapping_category, subcategory_1, subcategory_2):
    #remove mapping file samples that are not in the input BIOM table
    with open(mapping_fp, 'U') as map_f:
        md_map = MetadataMap.parseMetadataMap(map_f)
    md_map.filterSamples(loaded_biom.ids(axis='sample'), strict=True)

    if mapping_category not in md_map.CategoryNames:
        raise ValueError("category '%s' not found in mapping file "
                         "columns." % mapping_category)

    all_subcategories = md_map.getCategoryValues(md_map.sample_ids, mapping_category)

    if subcategory_1 not in all_subcategories:
        raise ValueError("subcategory_1 (-x) '%s' not found in selected "
                         "mapping file column." % subcategory_1)

    if subcategory_2 not in all_subcategories:
        raise ValueError("subcategory_2 (-y) '%s' not found in selected "
                         "mapping file column." % subcategory_2)

    if subcategory_2 == subcategory_1:
        raise ValueError("subcategory_1 (-x) must be different from subcategory_2 (-y)")

                    
def DA_fitZIG(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2):
   """perform metagenomeSeq's Zero Inflated Gaussian (ZIG) OTU differential abundance testing"""
   tmp_bt = load_table(input_path)
   tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
   check_mapping_file_category(tmp_bt, mapping_fp, mapping_category, subcategory_1, subcategory_2)
   tmp_bt.add_metadata(tmp_pmf, 'sample')

   with tempfile.NamedTemporaryFile(dir=get_qiime_temp_dir(),
                                    prefix='QIIME-differential-abundance-temp-table-',
                                    suffix='.biom') as temp_fh:
        temp_fh.write(tmp_bt.to_json('forR'))
        temp_fh.flush()
        run_fitZIG(temp_fh.name, out_path, mapping_category, subcategory_1, subcategory_2)

def multiple_file_DA_fitZIG(input_dir, output_dir, mapping_fp, mapping_category, subcategory_1, subcategory_2):
    """perform metagenomeSeq's Zero Inflated Gaussian (ZIG) OTU differential abundance test on a directory of raw abundance OTU matrices
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
        tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
        check_mapping_file_category(tmp_bt, mapping_fp, mapping_category, subcategory_1, subcategory_2)
        tmp_bt.add_metadata(tmp_pmf, 'sample')
        #make temporary json biom version - R currently does not have hdf5
        outfile = join(output_dir, 'fitZIG_DA_'+base_fname+'.txt')

        with tempfile.NamedTemporaryFile(dir=get_qiime_temp_dir(),
                                         prefix='QIIME-differential-abundance-temp-table-',
                                         suffix='.biom') as temp_fh:
            temp_fh.write(tmp_bt.to_json('forR'))
            temp_fh.flush()
            run_fitZIG(temp_fh.name, outfile, mapping_category, subcategory_1, subcategory_2) 


def run_fitZIG(input_path, out_path, mapping_category, subcategory_1, subcategory_2):
    """Run metagenomeSeq's fitZIG algorithm through Rscript
    """
    # set options
    command_args = ['-i %s -o %s -c %s -x %s -y %s' % (input_path, out_path, mapping_category, subcategory_1, subcategory_2)]
    # instantiate the object
    rsl = RExecutor(TmpDir=get_qiime_temp_dir())
    # run the app
    app_result = rsl(command_args=command_args, script_name='fitZIG.r')

    return app_result


def DA_DESeq2(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots):
    """perform DESeq2 negative binomial Wald differential abundance test on a raw abundance OTU matrix
    """
    tmp_bt = load_table(input_path)
    tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
    check_mapping_file_category(tmp_bt, mapping_fp, mapping_category, subcategory_1, subcategory_2)
    tmp_bt.add_metadata(tmp_pmf, 'sample')
    base_fname, ext = splitext(out_path)
    outfile_diagnostic = join(base_fname+'_diagnostic_plots.pdf') 

    with tempfile.NamedTemporaryFile(dir=get_qiime_temp_dir(),
                                     prefix='QIIME-differential-abundance-temp-table-',
                                     suffix='.biom') as temp_fh:
            temp_fh.write(tmp_bt.to_json('forR'))
            temp_fh.flush()
            run_DESeq2(temp_fh.name, out_path, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots, outfile_diagnostic) 


def multiple_file_DA_DESeq2(input_dir, output_dir, mapping_fp, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots):
    """perform DESeq2 negative binomial Wald differential abundance test on a directory of raw abundance OTU matrices
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
        tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
        check_mapping_file_category(tmp_bt, mapping_fp, mapping_category, subcategory_1, subcategory_2)
        tmp_bt.add_metadata(tmp_pmf, 'sample')
        outfile = join(output_dir, 'DESeq2_DA_'+base_fname+'.txt') 
        outfile_diagnostic = join(output_dir, 'DESeq2_diagnostic_plots_'+base_fname+'.pdf') 

        with tempfile.NamedTemporaryFile(dir=get_qiime_temp_dir(),
                                         prefix='QIIME-differential-abundance-temp-table-',
                                         suffix='.biom') as temp_fh:
            temp_fh.write(tmp_bt.to_json('forR'))
            temp_fh.flush()
            run_DESeq2(temp_fh.name, outfile, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots, outfile_diagnostic) 


def run_DESeq2(input_path, out_path, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots, outfile_diagnostic):
    """Run DESeq2 negative binomial Wald algorithm through Rscript
    """
    # set options
    if DESeq2_diagnostic_plots==True:
        command_args = ['-i %s -o %s -c %s -x %s -y %s -d %s -e %s' % (input_path, out_path, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots, outfile_diagnostic)]
    else:
        command_args = ['-i %s -o %s -c %s -x %s -y %s' % (input_path, out_path, mapping_category, subcategory_1, subcategory_2)]
    # instantiate the object
    rsl = RExecutor(TmpDir=get_qiime_temp_dir())
    # run the app
    app_result = rsl(command_args=command_args, script_name='DESeq2_nbinom.r')

    return app_result


def algorithm_list():
    """ returns list of differential abundance detection algorithms from qiime.differential_abundance
    """
    return ['metagenomeSeq_fitZIG', 'DESeq2_nbinom']
