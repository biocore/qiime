#!/usr/bin/env python
#file process_sff.py
from cogent.util.misc import app_path
from cogent.app.util import ApplicationNotFoundError
from cogent.parse.binary_sff import (
    parse_binary_sff, format_binary_sff, write_binary_sff, decode_accession,
    )
from os import listdir
from os.path import splitext, join, isfile, isdir, split
from cStringIO import StringIO
import itertools
import os
import subprocess

"""Converts directory of sff files into fasta and qual files.
"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight", "Greg Caporaso", "Jesse Stombaugh", "Kyle Bittinger"] 
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Release"


def _check_call(*args, **kwargs):
    """Run subprocess.check_call, sending stderr messages to /dev/null
    """
    kwargs['stderr'] = open(os.devnull, 'w')
    return subprocess.check_call(*args, **kwargs)


def _cumulative_sum(xs):
    """Return a list where each element is the sum up to that point.
    """
    cumsum = 0
    for x in xs:
        cumsum += x
        yield cumsum


def adjust_sff_cycles(sff_data, num_cycles):
    """Modify the number of cycles in a set of SFF data.

    This function reproduces the behavior of the -c option in Roche's
    sfffile program.
    """
    # TODO: Move to PyCogent
    num_flows = num_cycles * 4
    header, reads = sff_data

    h = header.copy()
    h['number_of_flows_per_read'] = num_flows
    h['header_length'] = num_flows + 40
    h['index_offset'] = 0
    h['index_length'] = 0
    h['flow_chars'] = 'TACG' * num_cycles

    read_clip_keys = [
        'clip_qual_left', 'clip_qual_right', 'clip_adapter_left',
        'clip_adapter_right',
        ]

    def adjust_read(read):
        r = read.copy()
        r['flowgram_values'] = read['flowgram_values'][:num_flows]
        enumerated_flow_indices = list(enumerate(
            _cumulative_sum(read['flow_index_per_base'])))

        # Brain teaser: find the largest base index having a flow
        # index less than num_flows
        num_bases = 6789
        for base_idx, flow_idx in reversed(enumerated_flow_indices):
            num_bases = base_idx + 1
            if flow_idx <= num_flows:
                break

        r['number_of_bases'] = num_bases
        r['flow_index_per_base'] = read['flow_index_per_base'][:num_bases]
        r['Bases'] = read['Bases'][:num_bases]
        r['quality_scores'] = read['quality_scores'][:num_bases]

        for key in read_clip_keys:
            if r[key] > num_bases:
                r[key] = num_bases
            
        return r

    return (h, itertools.imap(adjust_read, reads))


def format_binary_sff_as_fna(sff_file, output_file=None, qual=False):
    """Write a binary SFF file to an output file, in FASTA format.

    If no output file is provided, an in-memory file-like buffer is
    used (namely, a StringIO object).
    """
    # TODO: Move to PyCogent
    if output_file is None:
        output_file = StringIO()
    _, reads = parse_binary_sff(sff_file)
    for read in reads:
        output_file.write(format_read_as_fna(read, qual))
    return output_file


def format_read_as_fna(read, qual=False):
    """Format a single read from a binary SFF file as a FASTA record.

    If qual is True, output the qual values instead of the bases.
    """
    # TODO: Move to PyCogent
    out = StringIO()
    out.write('>%s' % read['Name'])

    # Roche uses 1-based indexing, where the right index is inclusive.
    # To transform to 0-based indices, where the right index is not
    # inclusive, we subtract 1 from the left index, but leave the
    # right index intact.

    start_idx = read['clip_qual_left'] - 1
    end_idx = read['clip_qual_right']

    # A surprising result is produced if the number of cycles are
    # adjusted such that no bases remain past clip_qual_left.  In the
    # clipping routine, the Roche software sets clip_qual_left to be
    # equal to the number of bases.  Using our indexing scheme, the
    # resulting sequence is of length 1 after clipping (one would
    # expect a length of 0).  We would fix this issue, if the effect
    # were not present in the output from Roche's sffinfo program.  We
    # retain this arguably incorrect behavior to be consistent with
    # the reference implementation.

    out.write(' length=%d' % (end_idx - start_idx))

    timestamp, _, region, location = decode_accession(read['Name'])
    out.write(' xy=%04d_%04d' % location)
    out.write(' region=%d' % region)
    out.write(' run=R_%d_%02d_%02d_%02d_%02d_%02d_' % timestamp)
    out.write('\n')

    if qual:
        scores = read['quality_scores'][start_idx:end_idx]
        out.write(' '.join(['%d' % s for s in scores]))
    else:
        bases = read['Bases'][start_idx:end_idx]
        out.write(bases)
    out.write('\n')
    return out.getvalue()


_MISSING_APP_MESSAGE = (
    "%s is not in $PATH. Is it installed? Have you added it to $PATH?")


def check_sffinfo():
    """Raise error if sffinfo is not in $PATH """
    if not app_path('sffinfo'):
        raise ApplicationNotFoundError(_MISSING_APP_MESSAGE % 'sffinfo')


def check_sfffile():
    """Raise error if sfffile is not in $PATH """
    if not app_path('sfffile'):
        raise ApplicationNotFoundError(_MISSING_APP_MESSAGE % 'sfffile')


def convert_Ti_to_FLX(sff_fp, output_fp, use_sfftools=False):
    """Converts Titanium SFF to FLX length reads."""
    if use_sfftools:
        check_sfffile()
        _check_call(
            ['sfffile', '-flx', '-o', output_fp, sff_fp],
            stdout=open(os.devnull, 'w'))
    else:
        header, reads = adjust_sff_cycles(parse_binary_sff(open(sff_fp), True), 100)
        write_binary_sff(open(output_fp, 'w'), header, reads)


def make_flow_txt(sff_fp, output_fp, use_sfftools=False):
    """Makes flowgram file from sff file."""
    if use_sfftools:
        check_sffinfo()
        _check_call(['sffinfo', sff_fp], stdout=open(output_fp, 'w'))
    else:
        format_binary_sff(open(sff_fp), open(output_fp, 'w'))


def make_fna(sff_fp, output_fp, use_sfftools=False):
    """Makes fna file from sff file."""
    if use_sfftools:
        check_sffinfo()
        _check_call(['sffinfo', '-s', sff_fp], stdout=open(output_fp, 'w'))
    else:
        format_binary_sff_as_fna(open(sff_fp), open(output_fp, 'w'))


def make_qual(sff_fp, output_fp, use_sfftools=False):
    """Makes qual file from sff file."""
    if use_sfftools:
        check_sffinfo()
        _check_call(['sffinfo', '-q', sff_fp], stdout=open(output_fp, 'w'))
    else:
        format_binary_sff_as_fna(open(sff_fp), open(output_fp, 'w'), qual=True)


def prep_sffs_in_dir(
    sff_dir, output_dir, make_flowgram=False, convert_to_flx=False,
    use_sfftools=False):
    """Converts all sffs in dir to fasta/qual.

    If convert_to_flx is True, each SFF file is first adjusted to
    contain flowgrams of length 100.

    If make_flowgram is true, the full text output is exported from
    the SFF file, in addition to the FASTA and qual files.

    If use_sfftools is True, the Roche programs sffinfo and sfffile
    will be called instead of the equivalent Python functions.
    """
    if isfile(sff_dir):
        # This is undocumented behavior, but we do support passing a
        # single file.  Do not check for sff extension; assume the
        # user knows what they are doing.
        filenames = [os.path.basename(sff_dir)]
        sff_dir = os.path.dirname(sff_dir)
    else:
        filenames = [x for x in os.listdir(sff_dir) if x.endswith('.sff')]

    for filename in filenames:
        sff_fp = join(sff_dir, filename)
        base_filename = splitext(filename)[0]
        base_output_fp = join(output_dir, base_filename)

        if convert_to_flx:
            base_output_fp = base_output_fp + '_FLX'
            sff_flx_fp = base_output_fp + '.sff'
            convert_Ti_to_FLX(sff_fp, sff_flx_fp, use_sfftools)
            # Converted sff file becomes basis for further processing
            sff_fp = sff_flx_fp

        make_fna(sff_fp, base_output_fp + '.fna', use_sfftools)
        make_qual(sff_fp, base_output_fp + '.qual', use_sfftools)
        if make_flowgram:
            make_flow_txt(sff_fp, base_output_fp + '.txt', use_sfftools)
