#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Daniel McDonald", "Greg Caporaso",
               "Justin Kuczynski", "Jens Reeder", "Catherine Lozupone",
               "Jai Ram Rideout", "Logan Knecht", "Michael Dwan",
               "Levi McCracken", "Damien Coy", "Yoshiki Vazquez Baeza",
               "Will Van Treuren", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


"""Contains general utility code in support of the Qiime project.

A lot of this might migrate into cogent at some point.
"""

import tempfile
import os
from os import getenv, listdir, close
from os.path import abspath, basename, exists, dirname, join, splitext, isfile
from collections import defaultdict
from gzip import open as gz_open
from sys import stderr
from copy import deepcopy
from datetime import datetime
from subprocess import Popen
from random import random
from itertools import repeat, izip
from functools import partial

from numpy import (array, zeros, shape, vstack, ndarray, asarray,
                   float, where, isnan, std, sqrt, ravel, mean, median,
                   sum as np_sum, nan, sort)
from numpy.ma import MaskedArray
from numpy.ma.extras import apply_along_axis

from biom.util import compute_counts_per_sample_stats, biom_open, HAVE_H5PY
from biom import load_table
from biom.exception import BiomException
from biom.table import Table

from cogent.parse.tree import DndParser
from cogent.cluster.procrustes import procrustes

from skbio.util import remove_files, create_dir
from skbio.format.sequences import format_fastq_record
from burrito.util import ApplicationError, CommandLineApplication, FilePath
from burrito.util import which
from skbio.sequence import DNASequence
from skbio.parse.sequences import parse_fasta
from skbio.parse.sequences import FastaIterator, FastqIterator

from bfillings.blast import Blastall, BlastResult
from bfillings.formatdb import (build_blast_db_from_fasta_path,
    build_blast_db_from_fasta_file)
from qiime_default_reference import (get_template_alignment,
    get_reference_sequences, get_reference_taxonomy,
    get_template_alignment_column_mask)


from qcli import make_option, qcli_system_call, parse_command_line_parameters

from qiime import __version__ as qiime_library_version
from qiime.parse import (parse_qiime_config_files,
                         parse_coords,
                         parse_newick,
                         fields_to_dict,
                         PhyloNode,
                         parse_mapping_file,
                         parse_denoiser_mapping,
                         parse_fastq,
                         mapping_file_to_dict)

# for backward compatibility - compute_seqs_per_library_stats has
# been removed in favor of biom.util.compute_counts_per_sample_stats,
# which has the same interface as the former
# qiime.util.compute_seqs_per_library_stats
compute_seqs_per_library_stats = compute_counts_per_sample_stats

# add a support message to script-raised errors
parse_command_line_parameters = partial(parse_command_line_parameters,
                                        error_suffix="\nIf you need help with "
                                                     "QIIME, see:\nhttp://help"
                                                     ".qiime.org\n")


class TreeMissingError(IOError):

    """Exception for missing tree file"""
    pass


class OtuMissingError(IOError):

    """Exception for missing OTU file"""
    pass


class AlignmentMissingError(IOError):

    """Exception for missing alignment file"""
    pass


class MissingFileError(IOError):
    pass


class FileFormatError(IOError):

    """Exception for wrong file format"""
    pass


class ScriptsDirError(IOError):
    """Exception for when the QIIME scripts directory cannot be found."""
    pass

class EmptyBIOMTableError(BiomException):
    """Exception for when an empty BIOM table is encountered."""
    pass

def make_safe_f(f, allowed_params):
    """Make version of f that ignores extra named params."""
    def inner(*args, **kwargs):
        if kwargs:
            new_kwargs = {}
            for k, v in kwargs.items():
                if k in allowed_params:
                    new_kwargs[k] = v
            return f(*args, **new_kwargs)
        return f(*args, **kwargs)
    return inner


class FunctionWithParams(object):

    """A FunctionWithParams is a replacement for the function factory.

    Specifically, the params that will be used in the __call__ method are
    available in a dict so you can keep track of them with the object
    itself.
    """
    Application = None
    Algorithm = None
    Citation = None
    Params = {}
    Name = 'FunctionWithParams'  # override in subclasses
    _tracked_properties = []  # properties tracked like params

    def __init__(self, params):
        """Return new FunctionWithParams object with specified params.

        Note: expect params to contain both generic and per-method (e.g. for
        cdhit) params, so leaving it as a dict rather than setting
        attributes.

        Some standard entries in params are:

        [fill in on a per-application basis]
        """
        self.Params.update(params)
        self._tracked_properties.extend(
            ['Application', 'Algorithm', 'Citation'])

    def __str__(self):
        """Returns formatted key-value pairs from params."""
        res = [self.Name + ' parameters:']
        for t in self._tracked_properties:
            res.append(t + ':' + str(getattr(self, t)))
        for k, v in sorted(self.Params.items()):
            res.append(str(k) + ':' + str(v))
        return '\n'.join(res)

    def writeLog(self, log_path):
        """Writes self.Params and other relevant info to supplied path."""
        f = open(log_path, 'w')
        f.write(str(self))
        f.close()

    def getResult(self, *args, **kwargs):
        """Gets result in __call__. Override in subclasses."""
        return None

    def formatResult(self, result):
        """Formats result as string (for whatever "result" means)."""
        return str(result)

    def writeResult(self, result_path, result):
        """Writes result to result_path. May need to format in subclasses."""
        f = open(result_path, 'w')
        f.write(self.formatResult(result))
        f.write('\n')
        f.close()

    def getTree(self, tree_source):
        """Returns parsed tree from putative tree source"""
        if isinstance(tree_source, PhyloNode):
            tree = tree_source  # accept tree object directly for tests
        elif tree_source:
            try:
                f = open(tree_source, 'U')
            except (TypeError, IOError):
                raise TreeMissingError(
                    "Couldn't read tree file at path: %s" %
                    tree_source)
            tree = parse_newick(f, PhyloNode)
            f.close()
        else:
            raise TreeMissingError(str(self.Name) +
                                   " is a phylogenetic metric, but no tree was supplied.")
        return tree

    def getData(self, data_source):
        """Returns data from putative source, which could be a path"""
        if isinstance(data_source, str):
            try:
                return eval(data_source)
            except (NameError, SyntaxError):
                try:
                    data_f = open(data_source, 'U')
                    data = data_f.read()
                    data_f.close()
                    try:
                        return eval(data)
                    except (NameError, SyntaxError, TypeError):
                        pass
                    return data
                except (IOError, NameError, TypeError):
                    pass
        # if we got here, either we didn't get a string or we couldn't read
        # the data source into any other kind of object
        return data_source

    def getBiomData(self, data):
        """returns a biom object regardless of whether path or object given"""
        try:
            if isfile(data):
                return load_table(data)
        except TypeError:
            if type(data) == Table:
                otu_table = data
                return otu_table
            else:
                raise TypeError('Data is neither a path to a biom table or a' +
                                ' biom table object.')

    def __call__(self, result_path=None, log_path=None,
                 *args, **kwargs):
        """Returns the result of calling the function using the params dict.

        Parameters:
        [fill in on a per-application basis]
        """
        result = self.getResult(*args, **kwargs)
        if log_path:
            self.writeLog(log_path)
        if result_path:
            self.writeResult(result_path, result)
        else:
            return result


def trim_fasta(fasta_lines, output_length):
    """trim fasta seqs to output_length bases """
    for seq_id, seq in parse_fasta(fasta_lines):
        yield '>%s\n%s\n' % (seq_id, seq[:output_length])


def get_qiime_project_dir():
    """ Returns the top-level QIIME directory
    """
    # Get the full path of util.py
    current_file_path = abspath(__file__)
    # Get the directory containing util.py
    current_dir_path = dirname(current_file_path)
    # Return the directory containing the directory containing util.py
    return dirname(current_dir_path)


def get_qiime_scripts_dir():
    """Return the directory containing QIIME scripts.

    The scripts directory is inferred from the location of
    print_qiime_config.py (equivalent to running ``which
    print_qiime_config.py``). Will raise a ``ScriptsDirError`` if the scripts
    directory cannot be determined.

    Note: this function will likely not work on Windows.

    """
    script_fp = which('print_qiime_config.py')

    if script_fp is None:
        raise ScriptsDirError("Could not find the directory containing QIIME "
                              "scripts. QIIME scripts must be accessible via "
                              "the PATH environment variable, and they must "
                              "be executable. Please ensure that you have a "
                              "valid QIIME installation (see the QIIME "
                              "Installation Guide: "
                              "http://qiime.org/install/install.html).")

    return dirname(script_fp)


def get_qiime_temp_dir():
    """ Returns the temp directory that should be used by QIIME scripts

    """
    return load_qiime_config()['temp_dir']


def load_qiime_config():
    """Return default parameters read in from file"""

    qiime_config_filepaths = []
    qiime_project_dir = get_qiime_project_dir()
    qiime_config_filepaths.append(
        qiime_project_dir + '/qiime/support_files/qiime_config')

    qiime_config_env_filepath = getenv('QIIME_CONFIG_FP')
    if qiime_config_env_filepath:
        qiime_config_filepaths.append(qiime_config_env_filepath)

    home_dir = getenv('HOME')
    if home_dir:
        qiime_config_home_filepath = home_dir + '/.qiime_config'
        qiime_config_filepaths.append(qiime_config_home_filepath)

    qiime_config_files = []
    for qiime_config_filepath in qiime_config_filepaths:
        if exists(qiime_config_filepath):
            qiime_config_files.append(open(qiime_config_filepath))

    qiime_config = parse_qiime_config_files(qiime_config_files)

    # For files that are defined in the qiime-default-reference package,
    # add values to the qiime_config if they haven't already been defined.
    qiime_config['pick_otus_reference_seqs_fp'] = \
        qiime_config['pick_otus_reference_seqs_fp'] or get_reference_sequences()

    qiime_config['pynast_template_alignment_fp'] = \
        qiime_config['pynast_template_alignment_fp'] or get_template_alignment()

    qiime_config['assign_taxonomy_reference_seqs_fp'] = \
        qiime_config['assign_taxonomy_reference_seqs_fp'] or get_reference_sequences()

    qiime_config['assign_taxonomy_id_to_taxonomy_fp'] = \
        qiime_config['assign_taxonomy_id_to_taxonomy_fp'] or get_reference_taxonomy()

    # Fall back to the system's temporary directory if one hasn't been defined.
    temp_dir = qiime_config['temp_dir'] or tempfile.gettempdir()

    # QIIME was historically written to assume that directory names end with a
    # path separator. Add one for safety.
    if not temp_dir.endswith(os.sep):
        temp_dir += os.sep
    qiime_config['temp_dir'] = temp_dir

    return qiime_config

def qiime_blast_seqs(seqs,
                     blast_constructor=Blastall,
                     blast_program='blastn',
                     blast_db=None,
                     refseqs=None,
                     refseqs_fp=None,
                     blast_mat_root=None,
                     params=None,
                     WorkingDir=None,
                     seqs_per_blast_run=1000,
                     is_protein=False,
                     HALT_EXEC=False):
    """Blast list of sequences.

    seqs: a list (or object with list-like interace) of (seq_id, seq)
     tuples (e.g., the output of parse_fasta)

    """

    assert blast_db or refseqs_fp or refseqs, \
        'Must provide either a blast_db or a fasta ' +\
        'filepath containing sequences to build one.'

    if refseqs_fp:
        blast_db, db_files_to_remove =\
            build_blast_db_from_fasta_path(refseqs_fp,
                                           output_dir=WorkingDir,
                                           is_protein=is_protein)
    elif refseqs:
        blast_db, db_files_to_remove =\
            build_blast_db_from_fasta_file(refseqs,
                                           output_dir=WorkingDir,
                                           is_protein=is_protein)
    else:
        db_files_to_remove = []

    if params is None:
        params = {}
    params["-d"] = blast_db
    params["-p"] = blast_program

    blast_app = blast_constructor(
        params=params,
        blast_mat_root=blast_mat_root,
        InputHandler='_input_as_seq_id_seq_pairs',
        WorkingDir=WorkingDir,
        SuppressStderr=True,
        HALT_EXEC=HALT_EXEC)

    current_seqs = []
    blast_results = BlastResult([])
    for seq in seqs:
        current_seqs.append(seq)
        if len(current_seqs) % seqs_per_blast_run == 0:
            if blast_results:
                blast_results.update(
                    BlastResult(blast_app(current_seqs)['StdOut']))
            else:
                blast_results = BlastResult(blast_app(current_seqs)['StdOut'])
            current_seqs = []

    # clean-up run: blast the remaining sequences
    blast_results.update(
        BlastResult(blast_app(current_seqs)['StdOut']))

    remove_files(db_files_to_remove)

    return blast_results


def qiime_blastx_seqs(seqs,
                      blast_constructor=Blastall,
                      blast_db=None,
                      refseqs=None,
                      refseqs_fp=None,
                      blast_mat_root=None,
                      params={},
                      WorkingDir=None,
                      seqs_per_blast_run=1000,
                      HALT_EXEC=False):
    """Blast list of sequences.

    seqs: a list (or object with list-like interace) of (seq_id, seq)
     tuples (e.g., the output of parse_fasta)

    """
    return qiime_blast_seqs(seqs,
                            blast_constructor=blast_constructor,
                            blast_program='blastx',
                            blast_db=blast_db,
                            refseqs=refseqs,
                            refseqs_fp=refseqs_fp,
                            blast_mat_root=blast_mat_root,
                            params={},
                            WorkingDir=WorkingDir,
                            seqs_per_blast_run=seqs_per_blast_run,
                            is_protein=True,
                            HALT_EXEC=HALT_EXEC)


def extract_seqs_by_sample_id(seqs, sample_ids, negate=False):
    """ Returns (seq id, seq) pairs if sample_id is in sample_ids """
    sample_ids = {}.fromkeys(sample_ids)

    if not negate:
        def f(s):
            return s in sample_ids
    else:
        def f(s):
            return s not in sample_ids

    for seq_id, seq in seqs:
        sample_id = seq_id.split('_')[0]
        if f(sample_id):
            yield seq_id, seq


def split_fasta_on_sample_ids(seqs):
    """ yields (sample_id, seq_id, seq) for each entry in seqs

        seqs: (seq_id,seq) pairs, as generated by parse_fasta

    """
    for seq_id, seq in seqs:
        yield (seq_id.split()[0].rsplit('_', 1)[0], seq_id, seq)
    return


def split_fasta_on_sample_ids_to_dict(seqs):
    """ return split_fasta_on_sample_ids as {sample_id: [(seq_id, seq), ], }

        seqs: (seq_id,seq) pairs, as generated by parse_fasta

    """
    result = {}
    for sample_id, seq_id, seq in split_fasta_on_sample_ids(seqs):
        try:
            result[sample_id].append((seq_id, seq))
        except KeyError:
            result[sample_id] = [(seq_id, seq)]
    return result


def write_seqs_to_fasta(fp, seqs, write_mode='w'):
    """Write seqs to fp with specified write mode ('a' or 'w')

        seqs: list of (seq_id,seq) tuples, as obtained from
         parse_fasta
    """
    with open(fp, write_mode) as f:
        for s in seqs:
            f.write('>%s\n%s\n' % (s))


def write_seqs_to_fastq(fp, seqs, write_mode='w'):
    """Write seqs to fp with specified write mode ('a' or 'w')

        seqs: list of (seq_id,seq,qual_id,qual) tuples
    """
    with open(fp, write_mode) as f:
        for s in seqs:
            f.write(format_fastq_record(s[0], s[1], s[3]))


def get_generated_by_for_biom_tables():
    """Returns a "generated by" string for use when creating BIOM tables

    Returns
    -------
    str
    """
    return "QIIME " + get_qiime_library_version()


def write_biom_table(biom_table, biom_table_fp, compress=True,
                     write_hdf5=HAVE_H5PY, table_type='OTU table'):
    """Writes a BIOM table to the specified filepath

    Parameters
    ----------
    biom_table : biom.Table
        The table object to write out
    biom_table_fp : str
        The path to the output file
    compress : bool, optional
        Defaults to ``True``. If True, built-in compression on the output HDF5
        file will be enabled. This option is only relevant if ``write_hdf5`` is
        ``True``.
    write_hdf5 : bool, optional
        Defaults to ``True`` if H5PY is installed and to ``False`` if H5PY is
        not installed. If ``True`` the output biom table will be written as an
        HDF5 binary file, otherwise it will be a JSON string.
    table_type : str, optional
        The Table.type value to set for the table before it is written. Note
        that this is a controlled vocabulary documented on biom-format.org.

    Raises
    ------
    EmptyBIOMTableError
        If ``biom_table.is_empty() == True``. 
    """
    if biom_table.is_empty():
        raise EmptyBIOMTableError(
            "Attempting to write an empty BIOM table to disk. "
            "QIIME doesn't support writing empty BIOM output files.")

    generated_by = get_generated_by_for_biom_tables()
    biom_table.type = table_type

    if write_hdf5:
        with biom_open(biom_table_fp, 'w') as biom_file:
            biom_table.to_hdf5(biom_file, generated_by, compress)
    else:
        with open(biom_table_fp, 'w') as biom_file:
            biom_table.to_json(generated_by, biom_file)


def split_sequence_file_on_sample_ids_to_files(seqs,
                                               file_type,
                                               output_dir,
                                               per_sample_buffer_size=500):
    """Split a demux sequence file into per-sample sequence file

        Parameters
        ----------
        seqs : file-like object
            The input file, which should be in the post-split_libraries format,
            where sequence identifiers are of the form sampleID_##
        file_type : {'fasta', 'fastq'}
            The type of the input file.
        output_dir : str
            Path to directory where output should be written, will be created
            if it doesn't exist
        per_sample_buffer_size : int, optional
            Defaults to 500. The number of sequences that will be accumulated
            for any given sample beore they are written out to disk.

        Notes
        -----
        - This function takes a buffered approach to writing files to avoid
          hitting errors arising from too many files being open when working
          with large numbers of samples ids (e.g. > 1024 on linux)
    """
    create_dir(output_dir)
    file_lookup = {}
    all_fps = []

    # Set these up here to reduce the number of checks that have to be
    # performed in the loop
    if file_type == 'fasta':
        iter_ = FastaIterator([seqs])
        sequence_writer = write_seqs_to_fasta
        data_getter = lambda x: (x[0], x[1])
    elif file_type == 'fastq':
        iter_ = FastqIterator([seqs])
        sequence_writer = write_seqs_to_fastq
        data_getter = lambda x: x
    else:
        raise ValueError("file_type must be either fasta or fastq")

    old_file_position = seqs.tell()
    seqs.seek(0)
    for record in iter_:
        seq_id = record['SequenceID']
        seq = record['Sequence']
        sample_id = seq_id.split()[0].rsplit('_', 1)[0]
        qual_id = record['QualID']
        qual = record['Qual']
        # grab or create the list corresponding to the current sample id
        try:
            current_seqs = file_lookup[sample_id][1]
        except KeyError:
            current_fp = join(output_dir, sample_id + '.%s' % file_type)
            all_fps.append(current_fp)
            if exists(current_fp):
                raise IOError(" %s already exists. Will not perform split "
                              "-- remove this file or specify a different "
                              "output directory." % current_fp)
            current_seqs = list()
            file_lookup[sample_id] = [current_fp, current_seqs]

        # append the current sequence to the current seqs list
        current_seqs.append(data_getter((seq_id, seq, qual_id, qual)))

        # compare current_seqs length to the buffer size, and write
        # if it has hit the buffer size
        if len(current_seqs) == per_sample_buffer_size:
            current_fp = file_lookup[sample_id][0]
            sequence_writer(current_fp, current_seqs, write_mode='a')
            # reset the current seqs buffer
            file_lookup[sample_id][1] = list()

    for current_fp, current_seqs in file_lookup.values():
        sequence_writer(current_fp, current_seqs, write_mode='a')

    seqs.seek(old_file_position)

    return all_fps


def median_absolute_deviation(x):
    """ compute the median of the absolute deviations from the median """
    x = asarray(x)
    median_x = median(x)
    median_abs_deviation = median(abs(x - median_x))
    return median_abs_deviation, median_x


def guess_even_sampling_depth(counts_per_sample, num_deviations=2.25):
    """ guess a depth for even sampling

        this is currently computed as the smallest seqs per sample
         count >= the median seqs per sample count - (2.25 * the median absolute
         deviation). 2.25 was chosen emprically by seeing how different values
         of num_deviations resulted in a choice that was similar to what
         I've chosen on several real OTU tables.
    """
    counts_per_sample.sort()
    median_abs_dev, median_count = \
        median_absolute_deviation(counts_per_sample)
    min_threshold = median_count - (num_deviations * median_abs_dev)
    for e in counts_per_sample:
        if e >= min_threshold:
            return e
    raise ValueError("No acceptable even sampling depth identified. " +
                     "It shouldn't be possible to get here, but just in case here's the " +
                     "counts per sample: %s" ' '.join(map(str, counts_per_sample)))


def raise_error_on_parallel_unavailable(qiime_config=None):
    """Raise error if no parallel QIIME bc user hasn't set jobs_to_start
    """
    if qiime_config is None:
        qiime_config = load_qiime_config()
    if 'jobs_to_start' not in qiime_config or \
       int(qiime_config['jobs_to_start']) < 2:
        raise RuntimeError("Parallel QIIME is not available. (Have you set" +
                           " jobs_to_start to greater than 1 in your qiime_config?")


def get_options_lookup():
    """ Return dict of commonly used options """
    qiime_config = load_qiime_config()
    result = {}
    result['fasta_as_primary_input'] =\
        make_option('-i', '--input_fasta_fp', type="existing_filepath",
                    help='path to the input fasta file')
    result['otu_table_as_primary_input'] =\
        make_option('-i', '--otu_table_fp', type="existing_filepath",
                    help='path to the input OTU table (i.e., the output from make_otu_table.py)')
    result['otu_map_as_primary_input'] =\
        make_option('-i', '--otu_map_fp', type="existing_filepath",
                    help='path to the input OTU map (i.e., the output from pick_otus.py)')
    result['log_fp'] =\
        make_option('-l', '--log_fp', type="new_filepath",
                    help='path to write the log file')
    result['input_fasta'] =\
        make_option('-f', '--input_fasta_fp', type="existing_filepath",
                    help='path to the input fasta file')
    result['output_dir'] =\
        make_option('-o', '--output_dir', type="new_dirpath",
                    help='path to the output directory')
    result['output_fp'] =\
        make_option('-o', '--output_fp', type="new_filepath",
                    help='the output filepath')
    result['output_biom_fp'] =\
        make_option('-o', '--output_biom_fp', type="new_filepath",
                    help='the output otu table in biom format (recommended extension: .biom)')
    result['mapping_fp'] =\
        make_option('-m', '--mapping_fp', type="existing_filepath",
                    help='the mapping filepath')

    # Define options used by the workflow scripts
    result['jobs_to_start_workflow'] =\
        make_option('-O', '--jobs_to_start', type='int',
                    help='Number of jobs to start. NOTE: you must also'
                    ' pass -a to run in parallel, this defines the number of'
                    ' jobs to be started if and only if -a is passed'
                    ' [default: %default]',
                    default=qiime_config['jobs_to_start'])

    # Define options used by the parallel scripts
    result['jobs_to_start'] =\
        make_option('-O', '--jobs_to_start', type='int',
                    help='Number of jobs to start [default: %default]',
                    default=qiime_config['jobs_to_start'])
    result['retain_temp_files'] =\
        make_option('-R', '--retain_temp_files', action='store_true',
                    help='retain temporary files after runs complete ' +
                    '(useful for debugging) [default: %default]',
                    default=False)
    result['suppress_submit_jobs'] =\
        make_option('-S', '--suppress_submit_jobs', action='store_true',
                    help='Only split input and write commands file - don\'t submit ' +
                    'jobs [default: %default]', default=False)
    result['poll_directly'] =\
        make_option('-T', '--poll_directly', action='store_true',
                    help='Poll directly for job completion rather than running ' +
                    'poller as a separate job. If -T is specified this script will ' +
                    'not return until all jobs have completed. [default: %default]',
                    default=False)
    result['cluster_jobs_fp'] =\
        make_option('-U', '--cluster_jobs_fp',
                    help='path to cluster jobs script (defined in qiime_config) ' +
                    ' [default: %default]',
                    default=qiime_config['cluster_jobs_fp'] or
                    'start_parallel_jobs.py')
    result['suppress_polling'] =\
        make_option('-W', '--suppress_polling', action='store_true',
                    help='suppress polling of jobs and merging of results ' +
                    'upon completion [default: %default]',
                    default=False)
    result['job_prefix'] =\
        make_option('-X', '--job_prefix', help='job prefix ' +
                    '[default: descriptive prefix + random chars]')
    result['seconds_to_sleep'] =\
        make_option('-Z', '--seconds_to_sleep', type='int',
                    help='Number of seconds to sleep between checks for run ' +
                    ' completion when polling runs [default: %default]',
                    default=qiime_config['seconds_to_sleep'] or 60)

    return result


def matrix_stats(headers_list, distmats):
    """does, mean, median, stdev on a series of (dis)similarity matrices

    takes a series of parsed matrices (list of headers, list of numpy 2d arrays)
    headers must are either row or colunm headers (those must be identical)
    outputs headers (list), means, medians, stdevs (all numpy 2d arrays)
    """

    if len(set(map(tuple, headers_list))) > 1:
        raise ValueError("error, not all input matrices have" +
                         " identical column/row headers")

    all_mats = array(distmats)  # 3d numpy array: mtx, row, col
    means = mean(all_mats, axis=0)
    medians = median(all_mats, axis=0)
    stdevs = std(all_mats, axis=0)

    return deepcopy(headers_list[0]), means, medians, stdevs


def convert_otu_table_relative(otu_table):
    """Convert the OTU table to relative abundances

    this method works on a parsed OTU table
    """
    sample_ids, otu_ids, otu_counts, consensus = otu_table
    otu_counts = asarray(otu_counts, float)
    otu_counts = otu_counts / otu_counts.sum()
    otu_counts = where(isnan(otu_counts), 0.0, otu_counts)
    return (sample_ids, otu_ids, otu_counts, consensus)


def convert_OTU_table_relative_abundance(otu_table):
    """convert the OTU table to have relative abundances rather than raw counts
    """
    output = []
    data_lines = []
    otu_ids = []
    tax_strings = []
    taxonomy = False
    for line in otu_table:
        line = line.strip().split('\t')
        if line[0].startswith('#OTU ID'):
            output.append('\t'.join(line))
            if line[-1] == 'Consensus Lineage':
                taxonomy = True
        elif line[0].startswith('#'):
            output.append('\t'.join(line))
        else:
            if taxonomy:
                vals = [float(i) for i in line[1:-1]]
                tax_strings.append(line[-1])
            else:
                vals = [float(i) for i in line[1:]]
                tax_string = None
            data = array(vals, dtype=float)
            data_lines.append(data)
            otu_ids.append(line[0])
    data_lines = array(data_lines)
    totals = sum(data_lines)
    new_values = []
    for i in data_lines:
        new_values.append(i / totals)
    for index, i in enumerate(new_values):
        line = [otu_ids[index]]
        line.extend([str(j) for j in i])
        if taxonomy:
            line.append(tax_strings[index])
        output.append('\t'.join(line))
    return output


def load_pcoa_files(pcoa_dir):
    """loads PCoA files from filepaths
    """
    support_pcoas = []
    pcoa_filenames = listdir(pcoa_dir)
    # ignore invisible files like .DS_Store
    pcoa_filenames = [fname for fname in pcoa_filenames if not
                      fname.startswith('.')]
    master_pcoa = open(join(pcoa_dir, pcoa_filenames[0]), 'U')
    master_pcoa = parse_coords(master_pcoa)
    for fname in pcoa_filenames:
        try:
            f = open(join(pcoa_dir, fname), 'U')
            pcoa_res = parse_coords(f)
            support_pcoas.append(pcoa_res)
            f.close()
        except IOError as err:
            stderr.write('error loading support pcoa ' + fname + '\n')
            exit(1)
    return master_pcoa, support_pcoas


def summarize_pcoas(master_pcoa, support_pcoas,
                    method='IQR', apply_procrustes=True):
    """returns the average PCoA vector values for the support pcoas

    Also returns the ranges as calculated with the specified method.
    The choices are:
        IQR: the Interquartile Range
        ideal fourths: Ideal fourths method as implemented in scipy
    """
    if apply_procrustes:
        # perform procrustes before averaging
        support_pcoas = [list(sp) for sp in support_pcoas]
        master_pcoa = list(master_pcoa)
        for i, pcoa in enumerate(support_pcoas):
            master_std, pcoa_std, m_squared = procrustes(
                master_pcoa[1], pcoa[1])
            support_pcoas[i][1] = pcoa_std
        master_pcoa[1] = master_std

    m_matrix = master_pcoa[1]
    m_eigvals = master_pcoa[2]
    m_names = master_pcoa[0]
    jn_flipped_matrices = []
    all_eigvals = []
    for rep in support_pcoas:
        matrix = rep[1]
        eigvals = rep[2]
        all_eigvals.append(eigvals)
        jn_flipped_matrices.append(_flip_vectors(matrix, m_matrix))
    matrix_average, matrix_low, matrix_high = _compute_jn_pcoa_avg_ranges(
        jn_flipped_matrices, method)
    # compute average eigvals
    all_eigvals_stack = vstack(all_eigvals)
    eigval_sum = np_sum(all_eigvals_stack, axis=0)
    eigval_average = eigval_sum / float(len(all_eigvals))
    return matrix_average, matrix_low, matrix_high, eigval_average, m_names


def _compute_jn_pcoa_avg_ranges(jn_flipped_matrices, method):
    """Computes PCoA average and ranges for jackknife plotting

    returns 1) an array of jn_averages
             2) an array of upper values of the ranges
            3) an array of lower values for the ranges

    method: the method by which to calculate the range
        IQR: Interquartile Range
        ideal fourths: Ideal fourths method as implemented in scipy
    """
    x, y = shape(jn_flipped_matrices[0])
    all_flat_matrices = [matrix.ravel() for matrix in jn_flipped_matrices]
    summary_matrix = vstack(all_flat_matrices)
    matrix_sum = np_sum(summary_matrix, axis=0)
    matrix_average = matrix_sum / float(len(jn_flipped_matrices))
    matrix_average = matrix_average.reshape(x, y)
    if method == 'IQR':
        result = matrix_IQR(summary_matrix)
        matrix_low = result[0].reshape(x, y)
        matrix_high = result[1].reshape(x, y)
    elif method == 'ideal_fourths':
        result = idealfourths(summary_matrix, axis=0)
        matrix_low = result[0].reshape(x, y)
        matrix_high = result[1].reshape(x, y)
    elif method == "sdev":
        # calculate std error for each sample in each dimension
        sdevs = zeros(shape=[x, y])
        for j in xrange(y):
            for i in xrange(x):
                vals = array([pcoa[i][j] for pcoa in jn_flipped_matrices])
                sdevs[i, j] = vals.std(ddof=1)
        matrix_low = -sdevs / 2
        matrix_high = sdevs / 2

    return matrix_average, matrix_low, matrix_high


def _flip_vectors(jn_matrix, m_matrix):
    """transforms PCA vectors so that signs are correct"""
    m_matrix_trans = m_matrix.transpose()
    jn_matrix_trans = jn_matrix.transpose()
    new_matrix = zeros(jn_matrix_trans.shape, float)
    for i, m_vector in enumerate(m_matrix_trans):
        jn_vector = jn_matrix_trans[i]
        disT = list(m_vector - jn_vector)
        disT = sum(map(abs, disT))
        jn_flip = jn_vector * [-1]
        disF = list(m_vector - jn_flip)
        disF = sum(map(abs, disF))
        if disT > disF:
            new_matrix[i] = jn_flip
        else:
            new_matrix[i] = jn_vector
    return new_matrix.transpose()


def IQR(x):
    """calculates the interquartile range of x

    x can be a list or an array

    returns min_val and  max_val of the IQR"""

    x.sort()
    # split values into lower and upper portions at the median
    odd = len(x) % 2
    midpoint = int(len(x) / 2)
    if odd:
        low_vals = x[:midpoint]
        high_vals = x[midpoint + 1:]
    else:  # if even
        low_vals = x[:midpoint]
        high_vals = x[midpoint:]
    # find the median of the low and high values
    min_val = median(low_vals)
    max_val = median(high_vals)
    return min_val, max_val


def matrix_IQR(x):
    """calculates the IQR for each column in an array
    """
    num_cols = x.shape[1]
    min_vals = zeros(num_cols)
    max_vals = zeros(num_cols)
    for i in range(x.shape[1]):
        col = x[:, i]
        min_vals[i], max_vals[i] = IQR(col)
    return min_vals, max_vals


def idealfourths(data, axis=None):
    """This function returns an estimate of the lower and upper quartiles of the data along
    the given axis, as computed with the ideal fourths. This function was taken
    from scipy.stats.mstat_extra.py (http://projects.scipy.org/scipy/browser/trunk/scipy/stats/mstats_extras.py?rev=6392)
    """
    def _idf(data):
        x = data.compressed()
        n = len(x)
        if n < 3:
            return [nan, nan]
        (j, h) = divmod(n / 4. + 5 / 12., 1)
        qlo = (1 - h) * x[j - 1] + h * x[j]
        k = n - j
        qup = (1 - h) * x[k] + h * x[k - 1]
        return [qlo, qup]
    data = sort(data, axis=axis).view(MaskedArray)
    if (axis is None):
        return _idf(data)
    else:
        return apply_along_axis(_idf, axis, data)


def isarray(a):
    """
    This function tests whether an object is an array
    """
    try:
        validity = isinstance(a, ndarray)
    except:
        validity = False

    return validity


def degap_fasta_aln(seqs):
    """degap a Fasta aligment.

    seqs: list of label,seq pairs
    """

    for (label, seq) in seqs:
        yield DNASequence(seq, id=label).degap()


def write_degapped_fasta_to_file(seqs, tmp_dir=get_qiime_temp_dir()):
    """ write degapped seqs to temp fasta file."""
    fd, tmp_filename = tempfile.mkstemp(dir=tmp_dir, prefix="degapped_",
                                        suffix=".fasta")
    close(fd)

    with open(tmp_filename, 'w') as fh:
        for seq in degap_fasta_aln(seqs):
            fh.write(seq.to_fasta())

    return tmp_filename


def get_diff_for_otu_maps(otu_map1, otu_map2):
    """return reads in two otu_maps that are not shared

    otu_map1, otu_map2: OTU to seqID mapping as dict of lists
    """

    otus1 = set(otu_map1.keys())
    otus2 = set(otu_map2.keys())
    ids1 = set([x for otu in otus1 for x in otu_map1[otu]])
    ids2 = set([x for otu in otus2 for x in otu_map2[otu]])

    return ids1 - ids2, ids2 - ids1


def compare_otu_maps(otu_map1, otu_map2, verbose=False):
    """compare two otu maps and compute fraction of

    otu_map1, otu_map2: OTU to seqID mapping as dict of lists
    """

    right = 0
    wrong = 0

    otus1 = set(otu_map1.keys())
    otus2 = set(otu_map2.keys())
    shared_otus = otus1.intersection(otus2)
    # check for equal members in shared OTUs
    for otu in shared_otus:
        members1 = set(otu_map1[otu])
        members2 = set(otu_map2[otu])

        right += len(members1 & members2)
        missing_in_2 = members1 - members2
        wrong += len(missing_in_2)
        if (verbose and len(missing_in_2) > 0):
            print "OTU id: %s" % otu
            print list(missing_in_2)
            print

    # process OTUs in 1 not in 2
    for otu in otus1 - shared_otus:
        wrong += len(otu_map1[otu])
        if verbose:
            print "OTU id: %s" % otu
            print list(otu_map1[otu])

    return float(wrong) / (right + wrong)


def compute_days_since_epoch(day, month, year):
    """ pass day, month, year to compute days since epoch (1/1/1970)

        Note that full years should always be provided: 09 is a
        different year than 2009!
    """
    d = datetime(int(year), int(month), int(day))
    epoch = datetime.utcfromtimestamp(0)
    return (d - epoch).days


def get_interesting_mapping_fields(mapping_data, mapping_headers):
    """ Returns headers for fields that are useful to color by in plots

        These fields are the ones that contain greater than one value
         and less values than the number of entries (so for example
         not SampleID)
    """
    result = []
    num_samples = len(mapping_data)
    num_cols = len(mapping_headers)
    transposed_data = array(mapping_data).T
    for header, datum in zip(mapping_headers, transposed_data):
        d = set(datum)
        len_d = len(d)
        if len_d > 1 and len_d < num_samples:
            result.append(header)
    return result


def flowgram_id_to_seq_id_map(seqs):
    """Map flowgram ids to sequence ids from a post-split_libraries fasta file
    """
    result = {}
    for id_, seq in seqs:
        fields = id_.split()
        seq_id = id_
        flowgram_id = fields[1]
        result[flowgram_id] = seq_id
    return result


def get_java_version():
    """Returns the Java version, or None if not installed"""
    o, e, exit_status = qiime_system_call("java -version")
    if exit_status != 0:
        return None

    # expect the first line of stderr to be 'java version "x.y.z_build"'
    e = e.strip().splitlines()
    version_line = e[0]
    if not version_line.startswith('java version'):
        return None
    else:
        return version_line.split()[-1].strip('"')

# retain qiime_system_call function name for backward compatibility
qiime_system_call = qcli_system_call


def get_qiime_library_version():
    """get QIIME version and the git SHA + current branch (if applicable)"""
    qiime_dir = get_qiime_project_dir()
    qiime_version = qiime_library_version

    # more information could be retrieved following this pattern
    sha_cmd = 'git --git-dir %s/.git rev-parse HEAD' % (qiime_dir)
    sha_o, sha_e, sha_r = qiime_system_call(sha_cmd)
    git_sha = sha_o.strip()

    branch_cmd = 'git --git-dir %s/.git rev-parse --abbrev-ref HEAD' %\
        (qiime_dir)
    branch_o, branch_e, branch_r = qiime_system_call(branch_cmd)
    git_branch = branch_o.strip()

    # validate the output from both command calls
    if is_valid_git_refname(git_branch) and is_valid_git_sha1(git_sha):
        return '%s, %s@%s' % (__version__, git_branch, git_sha[0:7])
    else:
        return '%s' % __version__


def is_valid_git_refname(refname):
    """check if a string is a valid branch-name/ref-name for git

    Input:
    refname: string to validate

    Output:
    True if 'refname' is a valid branch name in git. False if it fails to meet
    any of the criteria described in the man page for 'git check-ref-format',
    also see:

    http://www.kernel.org/pub/software/scm/git/docs/git-check-ref-format.html
    """
    if len(refname) == 0:
        return False

    # git imposes a few requirements to accept a string as a
    # refname/branch-name

    # They can include slash / for hierarchical (directory) grouping, but no
    # slash-separated component can begin with a dot . or end with the sequence
    # .lock
    if (len([True for element in refname.split('/')
            if element.startswith('.') or element.endswith('.lock')]) != 0):
        return False

    # They cannot have two consecutive dots .. anywhere
    if '..' in refname:
        return False

    # They cannot have ASCII control characters (i.e. bytes whose values are
    # lower than \040, or \177 DEL), space, tilde, caret ^, or colon : anywhere
    if len([True for refname_char in refname if ord(refname_char) < 40 or
            ord(refname_char) == 177]) != 0:
        return False
    if ' ' in refname or '~' in refname or '^' in refname or ':' in refname:
        return False

    # They cannot have question-mark ?, asterisk *, or open bracket [ anywhere
    if '?' in refname or '*' in refname or '[' in refname:
        return False

    # They cannot begin or end with a slash / or contain multiple consecutive
    # slashes
    if refname.startswith('/') or refname.endswith('/') or '//' in refname:
        return False

    # They cannot end with a dot ..
    if refname.endswith('.'):
        return False

    # They cannot contain a sequence @{
    if '@{' in refname:
        return False

    # They cannot contain a \
    if '\\' in refname:
        return False

    return True


def is_valid_git_sha1(hash):
    """check if a string is a valid git sha1 string

    Input:
    hash: string to validate

    Output:
    True if the string has 40 characters and is an hexadecimal number, False
    otherwise.

    """

    if len(hash) != 40:
        return False
    try:
        value = int(hash, 16)
    except ValueError:
        return False

    return True


def get_pynast_version():
    """Return PyNAST version string or None if PyNAST is not installed"""
    try:
        import pynast
        return pynast.__version__
    except ImportError:
        return None


def inflate_denoiser_output(
        centroid_seqs, singleton_seqs, denoiser_map, raw_seqs):
    """Expand denoiser fasta files based on denoiser map

        The inflation process works as follows: write each centroid
         sequence n times, where n is the number of reads in that
         cluster, and write each singleton once. While writing these
         out map back to original sequence identifiers.

        The seqs objects passed in are lists of (seq_id, seq) tuples,
         as returned from parse_fasta.


    """
    id_lookup = parse_denoiser_mapping(denoiser_map)
    flowgram_to_seq_id_lookup = flowgram_id_to_seq_id_map(raw_seqs)
    for id_, seq in centroid_seqs:
        # centroid headers look like
        #>FZTHQMS01E140G | cluster size: 4353
        id, cluster_size_str = id_.split(' | ')
        cluster_member_ids = id_lookup[id]
        for c in cluster_member_ids:
            yield flowgram_to_seq_id_lookup[c], seq

    for id_, seq in singleton_seqs:
        yield flowgram_to_seq_id_lookup[id_], seq

    return

# Functions for counting sequences in fasta files


def count_seqs(fasta_filepath, parser=parse_fasta):
    """ Count the sequences in fasta_filepath

        fasta_filepath: string indicating the full path to the file
    """
    # Open the file and pass it to py_count_seqs_from_file -- wrapping
    # this makes for easier unit testing
    return count_seqs_from_file(open(fasta_filepath, 'U'), parser=parser)


def count_seqs_from_file(fasta_file, parser=parse_fasta):
    """Return number of sequences in fasta_file (no format checking performed)

        fasta_file: an open file object

    """
    result = 0
    lens = []
    for record in parser(fasta_file):
        result += 1
        lens.append(len(record[1]))
    if result == 0:
        return result, None, None
    else:
        return result, mean(lens), std(lens)


def count_seqs_in_filepaths(fasta_filepaths, seq_counter=count_seqs):
    """ Wrapper to apply seq_counter to fasta_filepaths

        fasta_filepaths: list of one or more fasta filepaths
        seq_counter: a function which takes a single filepath
         and returns the count of the number of sequences
         (default: count_seqs) -- this is parameterized to
         facilitate unit testing
    """
    total = 0
    counts = []
    inaccessible_filepaths = []
    # iterate over the input files
    for fasta_filepath in fasta_filepaths:
        # if the file is actually fastq, use the fastq parser.
        # otherwise use the fasta parser
        if fasta_filepath.endswith('.fastq') or fasta_filepath.endswith('.fq'):
            parser = partial(parse_fastq, enforce_qual_range=False)
        elif fasta_filepath.endswith('.tre') or \
                fasta_filepath.endswith('.ph') or \
                fasta_filepath.endswith('.ntree'):
             # This is clunky, but really convenient bc
             # it lets us count tree tips with count_seqs.py
            def parser(f):
                t = DndParser(f, constructor=PhyloNode)
                return zip(t.iterTips(), repeat(''))
        else:
            parser = parse_fasta

        try:
            # get the count of sequences in the current file
            current_count = seq_counter(fasta_filepath, parser=parser)
            # store it
            counts.append((current_count, fasta_filepath))
            # and increment the total count
            total += current_count[0]
        except IOError:
            # if the file couldn't be open, keep track of the filepath
            inaccessible_filepaths.append(fasta_filepath)

    return counts, total, inaccessible_filepaths

# End functions for counting sequences in fasta files


def iseq_to_qseq_fields(line, barcode_in_header,
                        barcode_length, barcode_qual_c='b'):
    """ Split an Illumina sequence line into qseq fields"""
    record = line.strip().split(':')
    rec_0_1, rec_0_2 = record[0].split('_')
    rec_4_1, rec_4_23 = record[4].split('#')
    rec_4_2, rec_4_3 = rec_4_23.split('/')
    if barcode_in_header:
        barcode = rec_4_2[:barcode_length]
        sequence = record[5]
        barcode_qual = barcode_qual_c * barcode_length
        sequence_qual = record[6]
    else:
        barcode = record[5][:barcode_length]
        sequence = record[5][barcode_length:]
        barcode_qual = record[6][:barcode_length]
        sequence_qual = record[6][barcode_length:]
    return (rec_0_1, rec_0_2, record[1], record[2], record[3],
            rec_4_1, rec_4_2, rec_4_3), sequence, sequence_qual,\
        barcode, barcode_qual


def is_gzip(fp):
    """Checks the first two bytes of the file for the gzip magic number

    If the first two bytes of the file are 1f 8b (the "magic number" of a
    gzip file), return True; otherwise, return false.
    """
    return open(fp, 'rb').read(2) == '\x1f\x8b'


def gzip_open(fp):
    return gz_open(fp, 'rb')


def qiime_open(fp, permission='U'):
    """Wrapper to allow opening of gzipped or non-compressed files

    Read or write the contents of a file

    file_fp : file path
    permission : either 'r','w','a'

    If the file is binary, be sure to pass in a binary mode (append 'b' to
    the mode); opening a binary file in text mode (e.g., in default mode 'U')
    will have unpredictable results.
    """
    if is_gzip(fp):
        return gzip_open(fp)
    else:
        return open(fp, permission)


def make_compatible_distance_matrices(dm1, dm2, lookup=None):
    """ Intersect distance matrices and sort the values """
    dm1_ids = dm1[0]
    dm1_data = dm1[1]
    dm2_ids = dm2[0]
    dm2_data = dm2[1]

    if lookup:
        try:
            dm1_ids = [lookup[e] for e in dm1_ids]
            dm2_ids = [lookup[e] for e in dm2_ids]
        except KeyError as e:
            raise KeyError("All entries in both DMs must be in lookup if a "
                           "lookup is provided. Missing: %s" % str(e))
    order = [e for e in dm1_ids if e in dm2_ids]

    if len(order) == 0:
        return ([], []), ([], [])

    # store the intersected distance matrices here
    matrices = []

    # iterate over the distance matrices and identifiers to match the data
    # note that the order must be the same between the two matrices
    for ids, distance_matrix in [(dm1_ids, dm1_data), (dm2_ids, dm2_data)]:

        # the order is kept by getting the indices from this list
        indices = [ids.index(element) for element in order]

        # this matrix contains the matched up data
        out = distance_matrix[indices][:, indices]
        matrices.append(out)

    return (order, matrices[0]), (order, matrices[1])


def get_rdp_jarpath():
    """ Return jar file name for RDP classifier ($RDP_JAR_PATH)"""
    return getenv('RDP_JAR_PATH')


def expand_otu_ids(otu_map, otus_to_expand, ignore_missing=False):
    """From OTU map and otu ids, return seq ids represented by the OTUs
    """
    result = []
    for o in otus_to_expand:
        otu_id = o.split()[0]
        try:
            result += otu_map[otu_id]
        except KeyError:
            if ignore_missing:
                continue
            else:
                raise KeyError("OTU id not in OTU map: %s" % o.split()[0])
    return result
# This function (stderr) was pulled from the following website:
# http://www.java2s.com/Open-Source/Python/Math/SciPy/scipy/scipy/stats/stats.py.htm
# then modified to fit the purpose needed. Originally from Scipy.


def stderr(a, axis=0, ddof=1):
    """ Returns the estimated population standard error of the values in the
        passed array (i.e., N-1).  Axis can equal None (ravel array
        first), or an integer (the axis over which to operate).
    """
    a, axis = _chk_asarray(a, axis)
    return std(a, axis, ddof=1) / float(sqrt(a.shape[axis]))

# This function (_chk_asarray) was pulled from the following website:
# http://www.java2s.com/Open-Source/Python/Math/SciPy/scipy/scipy/stats/stats.py.htm
# then modified to fit the purpose needed. Originally from Scipy.


def _chk_asarray(a, axis):
    """ Converts a list into an numpy array """
    if axis is None:
        a = ravel(a)
        outaxis = 0
    else:
        a = asarray(a)
        outaxis = axis
    return a, outaxis


def subsample_fasta(input_fasta_fp,
                    output_fp,
                    percent_subsample):
    """ Writes random percent_sample of sequences from input fasta filepath

    input_fasta_fp: input fasta filepath
    output_fp: output fasta filepath
    percent_subsample: percent of sequences to write
    """

    input_fasta = open(input_fasta_fp, "U")

    output_fasta = open(output_fp, "w")

    for label, seq in parse_fasta(input_fasta):
        if random() < percent_subsample:
            output_fasta.write('>%s\n%s\n' % (label, seq))

    input_fasta.close()
    output_fasta.close()


def summarize_otu_sizes_from_otu_map(otu_map_f):
    """ Given an otu map file handle, summarizes the sizes of the OTUs

        This is useful for determining number of singletons, doubletons, etc
         from an OTU map.
    """
    result = {}
    for otu_id, seq_ids in fields_to_dict(otu_map_f).items():
        otu_size = len(seq_ids)
        try:
            result[otu_size] += 1
        except KeyError:
            result[otu_size] = 1

    result = sorted(result.items())
    return result


class MetadataMap():

    """This class represents a QIIME metadata mapping file.

    Public attributes:
        Comments - the comments associated with this metadata map (a list of
            strings)
    """

    req_header_prefix = ['SampleID', 'BarcodeSequence',
                         'LinkerPrimerSequence']
    req_header_suffix = ['Description']

    @staticmethod
    def parseMetadataMap(lines, case_insensitive=False):
        """Parses a QIIME metadata mapping file into a MetadataMap object.

        This static method is basically a factory that reads in the given
        metadata mapping file contents and returns a MetadataMap instance. This
        method is provided for convenience.

        Arguments:
            lines - a list of strings representing the file contents of a QIIME
                metadata mapping file
            case_insensitive - a boolean to uppercase non required mapping file
                headers
        """
        mapping_data, header, comments = parse_mapping_file(lines)
        if case_insensitive:
            req_header = set(MetadataMap.req_header_prefix + MetadataMap.req_header_suffix)
            for i in range(len(header)):
                if header[i] not in req_header:
                    header[i] = header[i].upper()
        mapping = mapping_file_to_dict(mapping_data, header)
        return MetadataMap(mapping, comments)

    @staticmethod
    def mergeMappingFiles(mapping_files, no_data_value='no_data',
                          case_insensitive=False):
        """ Merge list of mapping files into a single mapping file

            mapping_files: open file objects containing mapping data
            no_data_value: value to be used in cases where there is no
            mapping field associated with a sample ID (default: 'no_data')
        """
        func = partial(MetadataMap.parseMetadataMap,
                       case_insensitive=case_insensitive)
        metadata_maps = map(func, mapping_files)
        merged = sum(metadata_maps[1:], metadata_maps[0])
        merged.no_data_value = no_data_value
        return merged

    def __init__(self, sample_metadata, Comments):
        """Instantiates a MetadataMap object.

        Arguments:
            sample_metadata - the output of parse_mapping_file_to_dict(). It
                expects a python dict of dicts, where the top-level key is
                sample ID, and the inner dict maps category name to category
                value. This can be an empty dict altogether or the inner dict
                can be empty
            Comments - the output of parse_mapping_file_to_dict(). It expects a
                list of strings for the comments in the mapping file. Can be an
                empty list
        """
        self._metadata = sample_metadata
        self.Comments = Comments
        self.no_data_value = 'no_data'

    def __eq__(self, other):
        """Test this instance for equality with another.

        Note: This code was taken from http://stackoverflow.com/questions/
            390250/elegant-ways-to-support-equivalence-equality-in-python-
            classes.
        """
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        """Test this instance for inequality with another.

        Note: This code was taken from http://stackoverflow.com/questions/
            390250/elegant-ways-to-support-equivalence-equality-in-python-
            classes.
        """
        return not self.__eq__(other)

    def __str__(self):
        """Returns a tab-separated version of the mapping file

        Missing data will be filled in with the value set in the member
        variable no_data_value.

        Note: that required columns will be in the required positions, but
              optional columns may be re-ordered compared to the input data.
        """
        # this will hold the output lines while we are generating them
        output_lines = []

        # Build an ordered list of headers
        # 2. The optional columns in the mapping file
        headers_present = self._metadata.iteritems().next()[1].keys()
        optional_headers = list(set(headers_present) -
                                set(self.req_header_prefix +
                                    self.req_header_suffix))

        headers = (self.req_header_prefix + optional_headers +
                   self.req_header_suffix)

        output_lines.extend(self.Comments)
        output_lines.append('#' + '\t'.join(headers))

        for sample_id, data in self._metadata.iteritems():
            current_data = []

            # Get the first required columns
            current_data.append(sample_id)
            # skip the SampleID required header, since we get that from the
            # dict we are currently iterating over
            for header in self.req_header_prefix[1:]:
                if header not in data:
                    raise ValueError(
                        "Metadata mapping file is missing required column %r."
                        % header)
                current_data.append(data[header])

            # Get the optional columns; allow for None in these columns
            for header in optional_headers:
                value = self.no_data_value if data[header] is None else \
                    data[header]

                current_data.append(value)

            # get the last required columns
            for header in self.req_header_suffix:
                if header not in data:
                    raise ValueError(
                        "Metadata mapping file is missing required column %r."
                        % header)
                current_data.append(data[header])

            output_lines.append('\t'.join([str(x) for x in current_data]))

        return '\n'.join(output_lines) + '\n'

    def __add__(self, other):
        """Merges two mapping files

        Fills in None where there is no data (i.e., when one mapping file
        does not have a column in the other).  The comments from both mapping
        files will be concatenated.
        """
        # Make a defaultdict of defaultdicts, the latter of which returns
        # None when an key is not present
        merged_data = defaultdict(lambda: defaultdict(lambda: None))

        # We will keep track of all unique sample_ids and metadata headers
        # we have seen as we go
        all_sample_ids = set()
        all_headers = set()

        # add all values from self into the merged_data structure
        for sample_id, data in self._metadata.iteritems():
            all_sample_ids.add(sample_id)
            for header, value in data.iteritems():
                all_headers.add(header)
                merged_data[sample_id][header] = value

        # then add all data from other
        for sample_id, data in other._metadata.iteritems():
            all_sample_ids.add(sample_id)
            for header, value in data.iteritems():
                all_headers.add(header)
                # if the two mapping files have identical sample_ids and
                # metadata columns but have DIFFERENT values, raise a value
                # error
                if merged_data[sample_id][header] is not None and \
                        merged_data[sample_id][header] != value:
                    raise ValueError("Different values provided for %s for "
                                     "sample %s in different mapping files."
                                     % (header, sample_id))
                else:
                    merged_data[sample_id][header] = value

        # Now, convert what we have seen into a normal dict
        normal_dict = {}
        for sample_id in all_sample_ids:
            if sample_id not in normal_dict:
                normal_dict[sample_id] = {}

            for header in all_headers:
                normal_dict[sample_id][header] = \
                    merged_data[sample_id][header]

        # and create a MetadataMap object from it; concatenate comments
        return self.__class__(normal_dict, self.Comments + other.Comments)

    def getSampleMetadata(self, sample_id):
        """Returns the metadata associated with a particular sample.

        The metadata will be returned as a dict mapping category name to
        category value.

        Arguments:
            sample_id - the sample ID (string) to retrieve metadata for
        """
        return self._metadata[sample_id]

    def getCategoryValue(self, sample_id, category):
        """Returns the category value associated with a sample's category.

        The returned category value will be a string.

        Arguments:
            sample_id - the sample ID (string) to retrieve category information
                for
            category - the category name whose value will be returned
        """
        return self._metadata[sample_id][category]

    def getCategoryValues(self, sample_ids, category):
        """Returns all the values of a given category.

        The return categories will be a list.

        Arguments:
            sample_ids - An ordered list of sample IDs (i.e., from a distance
                matrix)
            category - the category name whose values will be returned
        """
        return [self._metadata[sid][category] for sid in sample_ids]

    def isNumericCategory(self, category):
        """Returns True if the category is numeric and False otherwise.

        A category is numeric if all values within the category can be
        converted to a float.

        Arguments:
            category - the category that will be checked
        """
        category_values = self.getCategoryValues(self.sample_ids, category)

        is_numeric = True
        for category_value in category_values:
            try:
                float(category_value)
            except ValueError:
                is_numeric = False
        return is_numeric

    def hasUniqueCategoryValues(self, category):
        """Returns True if the category's values are all unique.

        Arguments:
            category - the category that will be checked for uniqueness
        """
        category_values = self.getCategoryValues(self.sample_ids, category)

        is_unique = False
        if len(set(category_values)) == len(self.sample_ids):
            is_unique = True
        return is_unique

    def hasSingleCategoryValue(self, category):
        """Returns True if the category's values are all the same.

        For example, the category 'Treatment' only has values 'Control' for the
        entire column.

        Arguments:
            category - the category that will be checked
        """
        category_values = self.getCategoryValues(self.sample_ids, category)

        single_value = False
        if len(set(category_values)) == 1:
            single_value = True
        return single_value

    @property
    def sample_ids(self):
        """Returns the IDs of all samples in the metadata map.

        The sample IDs are returned as a list of strings in alphabetical order.
        """
        return sorted(self._metadata.keys())

    @property
    def CategoryNames(self):
        """Returns the names of all categories in the metadata map.

        The category names are returned as a list of strings in alphabetical
        order.
        """
        return sorted(self.getSampleMetadata(self.sample_ids[0]).keys()) \
            if len(self.sample_ids) > 0 else []

    def filterSamples(self, sample_ids_to_keep, strict=True):
        """Remove samples that are not in ``sample_ids_to_keep``.

        If ``strict=True``, a ``ValueError`` will be raised if any of the
        sample IDs in ``sample_ids_to_keep`` cannot be found in the metadata
        map.
        """
        for sid in self.sample_ids:
            if sid not in sample_ids_to_keep:
                del self._metadata[sid]

        if strict:
            extra_samples = set(sample_ids_to_keep) - set(self.sample_ids)

            if extra_samples:
                extra_samples_formatted = ', '.join(
                    map(lambda e: '"%s"' % e, extra_samples))
                raise ValueError(
                    "Could not find the following sample ID(s) in the "
                    "metadata mapping file: %s" % extra_samples_formatted)


class RExecutor(CommandLineApplication):

    """RExecutor application controller
       Runs R with a source script (from qiime/support_files/R)
    """
    _input_handler = '_input_as_path'
    _command = "R"
    _options = {}

    _R_parameters = {
        'flags': '--slave'
    }

    # The name of the R script (located under qiime/support_files/R/)
    _R_script = ''

    _parameters = {}
    _parameters.update(_options)
    _parameters.update(_R_parameters)

    def getHelp(self):
        """Returns documentation string"""
        help_str =\
            """
        Runs the specified r script using the specified command

        Outputs:
            The results of the r script that is ran
        """
        return help_str

    def __call__(self, command_args, script_name, verbose=False):
        """Run the specified r script using the commands_args

            returns a CommandLineAppResult object
        """
        self._R_script = script_name
        rscript = self._get_R_script_path()
        base_command = self._get_base_command()
        cd_command, base_command = base_command.split(';')
        cd_command += ';'
        R_source_dir = self._get_R_script_dir()

        # Build up the command, consisting of a BaseCommand followed by
        # input and output (file) specifications
        command = self._commandline_join(
            [cd_command, base_command,
                '--args',
                '--source_dir', R_source_dir,
             ] + command_args + [' < %s ' % (rscript)]
        )

        if self.HaltExec:
            raise AssertionError("Halted exec with command:\n%s" % command)

        # run command
        stdout, stderr, exit_status = qiime_system_call(command, shell=True)

        # Determine if error should be raised due to exit status of
        # appliciation
        if not self._accept_exit_status(exit_status):
            if exit_status == 2:
                raise ApplicationError(
                    'R library not installed:\nstdout:\n%s\nstderr:\n%s\n' %
                    (stdout, stderr))
            else:
                raise ApplicationError(
                    'Unacceptable application exit status: %d\ncommand: %s\n'
                    'stdout:\n%s\nstderr:\n%s\n' %
                    (exit_status, command, stdout, stderr))

        if verbose:
            print ('Command: %s\nstdout:\n%s\nstderr:\n%s\n' %
                   (command, stdout, stderr))

    # The methods below were taken from supervised_learning.py
    def _get_R_script_dir(self):
        """Returns the path to the qiime R source directory."""
        qiime_dir = get_qiime_project_dir()
        script_dir = join(qiime_dir, 'qiime', 'support_files', 'R')
        return script_dir

    def _get_R_script_path(self):
        """Returns the path to the R script to be executed."""
        return join(self._get_R_script_dir(), self._R_script)

    def _commandline_join(self, tokens):
        """Formats a list of tokens as a shell command."""
        commands = filter(None, map(str, tokens))
        return self._command_delimiter.join(commands).strip()

    def _accept_exit_status(self, exit_status):
        """ Return False to raise an error due to exit_status !=0."""
        if exit_status != 0:
            return False
        return True

    @property
    def RParameters(self):
        return self.__extract_parameters('R')

    def __extract_parameters(self, name):
        """Extracts parameters in self._<name>_parameters from self.Parameters.

        Allows the program to conveniently access a subset of user-
        adjusted parameters, which are stored in the Parameters
        attribute.

        Relies on the convention of providing dicts named according to
        "_<name>_parameters" and "_<name>_synonyms".  The main
        parameters object is expected to be initialized with the
        contents of these dicts.  This method will throw an exception
        if either convention is not adhered to.
        """
        parameters = getattr(self, '_' + name + '_parameters')
        result = Parameters(parameters)
        for key in result.keys():
            result[key] = self.Parameters[key]
        return result


def get_duplicates(fields):
    """ Returns duplicates out of a list

    Modified from stackoverflow.com example duplicate detection code
    http://stackoverflow.com/a/5420328

    fields:  list of elements to check for duplicates
    """
    cnt = {}
    for field in fields:
        try:
            cnt[field] += 1
        except KeyError:
            cnt[field] = 1
    return [key for key in cnt.keys() if cnt[key] > 1]


def duplicates_indices(fields):
    """ Gets dictionary of duplicates:locations in a list

    Modified from stackoverflow.com example duplicate detection code
    http://stackoverflow.com/a/5420328

    fields:  list of elements to check for duplicates
    """
    dup, ind = get_duplicates(fields), defaultdict(list)
    for i, v in enumerate(fields):
        if v in dup:
            ind[v].append(i)
    return ind


def head_gzip(fp, n=10):
    f = gzip_open(fp)
    for i in range(n):
        print f.readline(),


def add_filename_suffix(filepath, suffix):
    """Adds a suffix to the filepath, inserted before the file extension.

    Returns the new filepath string. For example, if filepath is 'foo.txt' and
    suffix is '_bar', 'foo_bar.txt' will be returned.

    Arguments:
        filepath - any filepath to append the suffix to (before the file
            extension, if it exists). Most useful if the filepath points to a
            file instead of a directory. The filepath isn't required to have
            an extension
    """
    root, extension = splitext(basename(filepath))
    return root + suffix + extension


def sync_biom_and_mf(pmf, bt):
    """Reduce mapping file dict and biom table to shared samples.

    Inputs:
     pmf - parsed mapping file from parse_mapping_file_to_dict (nested dict).
     bt - biom table object.
    Outputs are a bt and pmf that contain only shared samples and a set of
    samples that are not shared. If no samples are unshared this final output
    will be an empty set.
    """
    mf_samples = set(pmf)
    bt_samples = set(bt.ids())
    if mf_samples == bt_samples:
        # agreement, can continue without fear of breaking code
        return pmf, bt, set()
    else:
        shared_samples = mf_samples.intersection(bt_samples)
        # check that we shared something
        assert len(shared_samples) != 0, \
            "sync_biom_and_mf: No shared samples, no point in continuing."
        nonshared_samples = mf_samples.union(bt_samples) - shared_samples
        # remove samples that were in the mapping file but not biom file
        npmf = {k: v for k, v in pmf.items() if k in shared_samples}
        # remove samples in the biom table that were not in the mapping file

        def _f(sv, sid, smd):
            return sid in shared_samples
        nbt = bt.filter(_f, axis='sample')
    return npmf, nbt, nonshared_samples


def biom_taxonomy_formatter(bt, md_key):
    """Return md strings from bt using md_key in order of bt.observation_metadata

    There are multiple legacy formats for metadata encoding in biom formats
    including as lists, dicts, and strings. This function attempts to figure out
    what form the metadata is in and convert it into a single string. This
    function assumes that the metadata is encoded as a single format. It will
    break if some of the metadata is e.g., encoded as a dict and some as a list.

    Inputs:
     bt - biom table object
     md_key - string, the key to return the metadata from the biom table.
    Outputs a list of strings (in order of bt.observation_metadata entries) of
    metadata. If no metadata could be found using the given key the function
    will print a warning and return None.
    """
    if bt.metadata(axis='observation') is None:
        print 'Warning: No metadata in biom table. Won\'t alter calculations.'
        return None
    else:
        dtype = bt.metadata(axis='observation')[0][md_key]
    if isinstance(dtype, dict):
        data = []
        for md in bt.metadata(axis='observation'):
            tmp = []
            for k, v in md[md_key].iteritems():
                tmp.append('%s_%s' % (k, v))
            data.append(' '.join(tmp))
        # data = [' '.join(['%s_%s' % (k,v) for k,v in md[md_key].items()]) for \
        #     md in bt.observation_metadata]
        return map(str, data)
    elif isinstance(dtype, list):
        return (
            map(str, ['; '.join(md[md_key])
                      for md in bt.metadata(axis='observation')])
        )
    elif isinstance(dtype, (str, unicode)):
        return map(str, [md[md_key] for md in bt.metadata(axis='observation')])
    else:
        print ('Metadata format could not be determined or metadata key (%s) ' +
               'was incorrect. Metadata will not be returned.') % md_key
        return None


def invert_dict(d):
    """Returns inverse of d, setting keys to values and values to list of keys.

    Note that each value will _always_ be a list, even if one item.

    Can be invoked with anything that can be an argument for dict(), including
    an existing dict or a list of tuples. However, keys are always appended in
    arbitrary order, not the input order.

    WARNING: will fail if any values are unhashable, e.g. if they are dicts or
    lists.

    Ported from PyCogent's cogent.util.misc.InverseDictMulti.
    """
    if isinstance(d, dict):
        temp = d
    else:
        temp = dict(d)
    result = {}
    for key, val in temp.iteritems():
        if val not in result:
            result[val] = []
        result[val].append(key)
    return result
