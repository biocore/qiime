#!/usr/bin/env python
# File created on 24 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from numpy import array, in1d
from itertools import product

from skbio.parse.sequences import parse_fasta
from skbio.util import create_dir

from qiime.parse import parse_mapping_file
from qiime.format import format_mapping_file


def make_field_value_list(headers, field, mdata):
    '''Return sorted list of unique values field takes in mdata.

    Parameters
    ----------
    headers : list
        Strings that are the header fields in a mapping file. Usually derived
        from parse_mapping_file.
    field : str
        Header of interest in headers.
    mdata : np.array
        2-d array, containing data from mapping file cast as array of strings.
        Usually derived from parse_mapping_file.

    Returns
    -------
    list
        Sorted list of unique values found in mapping field.

    Notes
    -----
    Function returns a sorted list rather than set because it keeps the order in
    memory the same allowing test code to work more easily. Performance cost is 
    tiny.

    Examples
    --------
    >>> from qiime.split import make_field_value_list
    >>> from numpy import array
    >>> headers = ['color', 'temp', 'size']
    >>> mdata = array([['s0', 'blue', 'hot', '13'],
                       ['s1', 'blue', 'cold', '1'],
                       ['s2', 'green', 'cold', '12'],
                       ['s3', 'cyan', 'hot', '1'],
                       ['s4', 'blue', '0', '0']], 
                      dtype='|S5')
    >>> make_field_value_list(headers, 'color', mdata)
    ['blue', 'cyan', 'green']
    '''
    return sorted(set(mdata[:, headers.index(field)]))

def make_field_set_iterable(headers, fields, mdata):
    '''Return product of lists of unique values in order of the passed fields.

    Parameters
    ----------
    headers : list
        Strings that are the header fields in a mapping file. Usually derived
        from parse_mapping_file.
    fields : list
        List of strings, headers of interest in headers.
    mdata : np.array
        2-d array, containing data from mapping file cast as array of strings.
        Usually derived from parse_mapping_file.

    Returns
    -------
    generator
        Generator that yields successive elements of the Cartesian product of
        the input lists.

    Examples
    --------
    >>> from qiime.split import make_field_set_iterable
    >>> from numpy import array
    >>> headers = ['color', 'temp', 'size']
    >>> mdata = array([['s0', 'blue', 'hot', '13'],
                       ['s1', 'blue', 'cold', '1'],
                       ['s2', 'green', 'cold', '12'],
                       ['s3', 'cyan', 'hot', '1'],
                       ['s4', 'blue', '0', '0']], 
                      dtype='|S5')
    >>> list(make_field_set_iterable(['color', 'temp'], headers, mdata)
    [('blue', '0'),
     ('blue', 'cold'),
     ('blue', 'hot'),
     ('cyan', '0'),
     ('cyan', 'cold'),
     ('cyan', 'hot'),
     ('green', '0'),
     ('green', 'cold'),
     ('green', 'hot')]
    '''
    return product(*[make_field_value_list(headers, f, mdata) for f in fields])

def make_non_empty_sample_lists(fields, headers, mdata):
    '''Return non-empty sample lists for corresponding field value sets.

    Parameters
    ----------
    headers : list
        Strings that are the header fields in a mapping file. Usually derived
        from parse_mapping_file.
    fields : list
        List of strings, headers of interest in headers.
    mdata : np.array
        2-d array, containing data from mapping file cast as array of strings.
        Usually derived from parse_mapping_file.

    Returns
    -------
    sample_groups : list
        A list of arrays where each array contains the samples that had fields
        equal to the given values in value_groups. Empty arrays are not
        returned.
    value_groups : list
        A list of tuples representing the values that the fields of interest
        took for the corresponding sample group in sample_groups.

    Examples
    --------
    >>> from qiime.split import make_field_set_iterable
    >>> from numpy import array
    >>> headers = ['color', 'temp', 'size']
    >>> mdata = array([['s0', 'blue', 'hot', '13'],
                       ['s1', 'blue', 'cold', '1'],
                       ['s2', 'green', 'cold', '12'],
                       ['s3', 'cyan', 'hot', '1'],
                       ['s4', 'blue', '0', '0']], 
                      dtype='|S5')
    >>> sgs, vgs = make_sample_lists(['color', 'temp'], headers, mdata)
    >>> sgs
    [array(['s4'], dtype='|S5'),
     array(['s1'], dtype='|S5'),
     array(['s0'], dtype='|S5'),
     array(['s3'], dtype='|S5'),
     array(['s2'], dtype='|S5')]
    >>> # notice that since there were no combinations of a sample that was 
    >>> # both cyan and cold, it is not included in the output.
    >>> vgs
    [('blue', '0'),
     ('blue', 'cold'),
     ('blue', 'hot'),
     ('cyan', 'hot'),
     ('green', 'cold')]
    '''
    fsi = make_field_set_iterable(headers, fields, mdata)

    # subset the data columns so we can operate on a smaller array. metadata
    # is an array of just the fields of the mapping file data that we are 
    # interested in.
    samples = mdata[:, 0]
    f_inds = [headers.index(i) for i in fields]
    metadata = mdata[:, f_inds]
    
    sample_groups = []
    value_groups = []
    for value_set in fsi:
        rows, = (metadata == value_set).all(1).nonzero()
        if rows.size > 0:
            sample_groups.append(samples[rows])
            value_groups.append(value_set)
        else:
            pass

    return sample_groups, value_groups

def subset_mapping_data(mdata, samples_of_interest):
    '''Remove rows of mdata that are not from samples_of_interest.

    Parameters
    ----------
    mdata : np.array
        2-d array, containing data from mapping file cast as array of strings.
        Usually derived from parse_mapping_file.
    samples_of_interest : list
        A list of strings that are a strict subset of the samples found in the
        first column of mdata.

    Returns
    -------
    subset of mdata
    
    Examples
    --------
    >>> from qiime.split import subset_mapping_data
    >>> from numpy import array
    >>> mdata = array([['s0', 'blue', 'hot', '13'],
                       ['s1', 'blue', 'cold', '1'],
                       ['s2', 'green', 'cold', '12'],
                       ['s3', 'cyan', 'hot', '1'],
                       ['s4', 'blue', '0', '0']], 
                      dtype='|S5')
    >>> subset_mapping_data(mdata, ['s0', 's2'])
    array([['s0', 'blue', 'hot', '13'],
           ['s2', 'green', 'cold', '12']], 
          dtype='|S5')
    '''
    return mdata[in1d(mdata[:, 0], samples_of_interest)]


def split_fasta(infile, seqs_per_file, outfile_prefix, working_dir=''):
    """ Split infile into files with seqs_per_file sequences in each

        infile: list of fasta lines or open file object
        seqs_per_file: the number of sequences to include in each file
        out_fileprefix: string used to create output filepath - output filepaths
         are <out_prefix>.<i>.fasta where i runs from 0 to number of output files
        working_dir: directory to prepend to temp filepaths (defaults to
         empty string -- files written to cwd)

        List of output filepaths is returned.

    """
    if seqs_per_file <= 0:
        raise ValueError("seqs_per_file must be > 0!")

    seq_counter = 0
    out_files = []
    if working_dir and not working_dir.endswith('/'):
        working_dir += '/'
        create_dir(working_dir)

    for seq_id, seq in parse_fasta(infile):
        if seq_counter == 0:
            current_out_fp = '%s%s.%d.fasta' \
                % (working_dir, outfile_prefix, len(out_files))
            current_out_file = open(current_out_fp, 'w')
            out_files.append(current_out_fp)
        current_out_file.write('>%s\n%s\n' % (seq_id, seq))
        seq_counter += 1

        if seq_counter == seqs_per_file:
            current_out_file.close()
            seq_counter = 0

    if not current_out_file.closed:
        current_out_file.close()

    return out_files
