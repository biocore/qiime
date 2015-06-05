#!/usr/bin/env python

"""Preprocess 454 sequencing data."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from itertools import imap
from os import remove, close
from random import sample
from collections import defaultdict
from string import lowercase
from tempfile import mkstemp

from skbio.parse.sequences import parse_fasta
from bfillings.denoiser import (Flowgram, build_averaged_flowgram,
                             lazy_parse_sff_handle, build_prefix_map)

from qiime.util import load_qiime_config
from qiime.denoiser.cluster_utils import submit_jobs
from qiime.denoiser.flowgram_filter import cleanup_sff,\
    truncate_flowgrams_in_SFF, extract_barcodes_from_mapping
from qiime.denoiser.utils import squeeze_seq, make_stats, get_representatives,\
    wait_for_file, store_mapping, invert_mapping, cat_sff_files, files_exist,\
    read_denoiser_mapping, get_denoiser_data_dir, write_sff_header

STANDARD_BACTERIAL_PRIMER = "CATGCTGCCTCCCGTAGGAGT"


def make_tmp_name(length=8):
    """Returns a random string of specified length.

    length: length of random string
    """
    return ("".join(sample(list(lowercase), length)))


def sample_mapped_keys(mapping, min_coverage=50):
    """sample up to min_coverage keys for each key in mapping.

    mapping: dictionary of lists.

    Note: key is always included in sample
    """
    if min_coverage == 0:
        return {}
    sample_keys = {}
    for key in mapping.keys():
        if (min_coverage > 1):
            sample_keys[key] = sample(mapping[key],
                                      min(min_coverage - 1, len(mapping[key])))
        else:
            sample_keys[key] = []
        sample_keys[key].append(key)  # always include the centroid
    return sample_keys


def build_averaged_flowgrams(mapping, sff_fp,
                             min_coverage=50, out_fp=None):
    """Build averaged flowgrams for each cluster in mapping.

    mapping: a cluster mapping as dictionary of lists

    sff_fp: pointer to sff.txt file, must be consistent with  mapping

    min_coverage: number of flowgrams to average over for each cluster

    out_fp: ouput file name

    NOTE: This function has no test code, since it is mostly IO around tested functions
    """

    l = len(mapping)
    (flowgrams, header) = lazy_parse_sff_handle(open(sff_fp))
    # update some values in the sff header
    header["# of Reads"] = l
    header["Index Length"] = "NA"

    if (out_fp):
        out_filename = out_fp
    else:
        fd, out_filename = mkstemp(dir="/tmp/",
                                  prefix="prefix_dereplicated",
                                  suffix=".sff.txt")
        close(fd)
    outhandle = open(out_filename, "w")

    # write out reduced flogram set
    write_sff_header(header, outhandle)

    seqs = {}
    # get a random sample for each cluster
    sample_keys = sample_mapped_keys(mapping, min_coverage)
    for ave_f, id in _average_flowgrams(mapping, flowgrams, sample_keys):
        outhandle.write(ave_f.createFlowHeader() + "\n")
        ave_f.Bases = ave_f.toSeq()
        seqs[id] = ave_f.Bases

    outhandle.close()
    return(out_filename, seqs)


def _average_flowgrams(mapping, flowgrams, sample_keys):
    """average flowgrams according to cluster mapping.

    mapping: a dictionary of lists as cluster mapping

    flowgrams:  an iterable flowgram source, all flowgram ids from this source must be in the mapping

    sample_keys: the keys that should be averaged over for each cluster.
    """

    # accumulates flowgram for each key until sample for this key is empty
    flows = defaultdict(list)
    invert_map = invert_mapping(mapping)
    for f in flowgrams:
        key = invert_map[f.Name]
        samples = sample_keys[key]
        if (f.Name in samples):
            flows[key].append(f.flowgram)
            samples.remove(f.Name)
            if (len(samples) == 0):
                # we gathered all sampled flowgrams for this cluster,
                # now average
                ave_flowgram = build_averaged_flowgram(flows[key])
                ave_f = Flowgram(ave_flowgram, Name=key)

                del(flows[key])
                yield ave_f, key


def prefix_filter_flowgrams(flowgrams, squeeze=False):
    """Filters flowgrams by common prefixes.

    flowgrams: iterable source of flowgrams

    squeeze: if True, collapse all poly-X to X

    Returns prefix mapping.
    """

    # collect flowgram sequences
    if squeeze:
        seqs = imap(
            lambda f: (f.Name, squeeze_seq(str(f.toSeq(truncate=True)))),
            flowgrams)
    else:
        seqs = imap(lambda f: (f.Name, str(f.toSeq(truncate=True))), flowgrams)
    # equivalent but more efficient than
    #seqs = [(f.Name, str(f.toSeq(truncate=True))) for f in flowgrams]

    # get prefix mappings
    mapping = build_prefix_map(seqs)
    l = len(mapping)
    orig_l = sum([len(a) for a in mapping.values()]) + l

    return (l, orig_l, mapping)


def print_rep_seqs(mapping, seqs, out_fp):
    """Print the cluster seeds of a mapping to out_fp.

    mapping: a cluster mapping

    seqs: a list of seqs contained in the mapping

    out_fp: output directory
    """
    out_fh = open(out_fp + "/prefix_dereplicated.fasta", "w")
    for s in (get_representatives(mapping, seqs.iteritems())):
        out_fh.write(s.to_fasta())
    out_fh.close()


def preprocess(sff_fps, log_fh, fasta_fp=None, out_fp="/tmp/",
               verbose=False, squeeze=False,
               primer=STANDARD_BACTERIAL_PRIMER):
    """Quality filtering and truncation of flowgrams, followed by denoiser phase I.

    sff_fps: List of paths to flowgram files

    log_fh: log messages are written to log_fh if it is set to something else than None

    fasta_fp: Path to fasta file, formatted as from split_libraries.py.
              This files is used to filter the flowgrams in sff_fps. Only reads in
              fasta_fp are pulled from sff_fps.

    out_fp: path to output directory

    verbose: a binary verbose flag

    squeeze: a flag that controls if sequences are squeezed before phase I.
             Squeezing means consecutive identical nucs are collapsed to one.

    primer: The primer sequences of the amplification process. This seq will be
            removed from all reads during the preprocessing
    """
    flowgrams, header = cat_sff_files(map(open, sff_fps))

    if(fasta_fp):
        # remove barcodes and sequences tossed by split_libraries, i.e. not in
        # fasta_fp
        labels = imap(lambda a_b: a_b[0], parse_fasta(open(fasta_fp)))
        barcode_mapping = extract_barcodes_from_mapping(labels)
        (trunc_sff_fp, l) = truncate_flowgrams_in_SFF(flowgrams, header,
                                                      outdir=out_fp,
                                                      barcode_mapping=barcode_mapping,
                                                      primer=primer)
        if verbose:
            log_fh.write(
                "Sequences in barcode mapping: %d\n" %
                len(barcode_mapping))
            log_fh.write("Truncated flowgrams written: %d\n" % l)
    else:
        # just do a simple clean and truncate
        (clean_sff_fp, l) = cleanup_sff(flowgrams, header, outdir=out_fp)
        if verbose:
            log_fh.write("Cleaned flowgrams written: %d\n" % l)
        flowgrams, header = lazy_parse_sff_handle(open(clean_sff_fp))
        (trunc_sff_fp, l) = truncate_flowgrams_in_SFF(flowgrams, header,
                                                      outdir=out_fp, primer=primer)
        if verbose:
            log_fh.write("Truncated flowgrams written: %d\n" % l)
        remove(clean_sff_fp)

    if (l == 0):
        raise ValueError("No flowgrams left after preprocesing.\n" +
                         "Check your primer sequence")

    # Phase I - cluster seqs which are exact prefixe
    if verbose:
        log_fh.write("Filter flowgrams by prefix matching\n")

    (flowgrams, header) = lazy_parse_sff_handle(open(trunc_sff_fp))
    l, orig_l, mapping =\
        prefix_filter_flowgrams(flowgrams, squeeze=squeeze)

    averaged_sff_fp, seqs = build_averaged_flowgrams(mapping, trunc_sff_fp,
                                                     min_coverage=1,
                                                     # averaging produces too good flowgrams
                                                     # such that the greedy clustering clusters too much.
                                                     # Use the cluster centroid
                                                     # instead by using
                                                     # min_coverage 1
                                                     out_fp=out_fp + "/prefix_dereplicated.sff.txt")
    remove(trunc_sff_fp)
    if verbose:
        log_fh.write("Prefix matching: removed %d out of %d seqs\n"
                     % (orig_l - l, orig_l))
        log_fh.write("Remaining number of sequences: %d\n" % l)
        log_fh.write(make_stats(mapping) + "\n")

    # print representative sequences and mapping
    print_rep_seqs(mapping, seqs, out_fp)
    store_mapping(mapping, out_fp, "prefix")
    return (averaged_sff_fp, l, mapping, seqs)


def preprocess_on_cluster(sff_fps, log_fp, fasta_fp=None, out_fp="/tmp/",
                          squeeze=False, verbose=False,
                          primer=STANDARD_BACTERIAL_PRIMER):
    """Call preprocess via cluster_jobs_script on the cluster.

    sff_fps: List of paths to flowgram files.

    log_fp: path to log file

    fasta_fp: Path to fasta file, formatted as from split_libraries.py.
              This files is used to filter the flowgrams in sff_fps. Only reads in
              fasta_fp are pulled from sff_fps.

    out_fp: path to output directory

    verbose: a binary verbose flag

    squeeze: a flag that controls if sequences are squeezed before phase I.
             Squeezing means consecutive identical nucs are collapsed to one.

    primer: The primer sequences of the amplification process. This seq will be
            removed from all reads during the preprocessing
    """
    cmd = "denoiser_preprocess.py -i %s -l %s -o %s" %\
        (",".join(sff_fps), log_fp, out_fp)
    if (fasta_fp):
        cmd += " -f %s" % fasta_fp
    if(squeeze):
        cmd += " -s"
    if verbose:
        cmd += " -v"
    if primer:
        cmd += " -p %s" % primer

    submit_jobs([cmd], "pp_" + make_tmp_name(6))

    wait_for_file(out_fp + "/prefix_mapping.txt", 10)


def read_preprocessed_data(out_fp="/tmp/"):
    """Read data of a previous preprocessing run.

    out_fp: output directory of previous preprocess run.
            Supposed to contain two files:
              - prefix_dereplicated.fasta
              - prefix_mapping.txt
    """
    # read mapping, and extract seqs
    # mapping has fasta_header like this:
    #  > id:   count

    seqs = dict([(a.split(':')[0], b) for (a, b) in
                (parse_fasta(open(out_fp + "/prefix_dereplicated.fasta")))])

    mapping = read_denoiser_mapping(open(out_fp + "/prefix_mapping.txt"))

    return(out_fp + "/prefix_dereplicated.sff.txt", len(mapping), mapping, seqs)
