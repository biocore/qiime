#!/usr/bin/env python
from __future__ import division

"""Various helper functions."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

import sys
from os import remove, makedirs, access, X_OK, R_OK, close
from os.path import exists, isdir
from collections import defaultdict
from re import sub
from time import sleep
from socket import error
from itertools import chain
from subprocess import Popen, PIPE, STDOUT
import pickle
from tempfile import mkstemp

from skbio.sequence import BiologicalSequence
from burrito.util import ApplicationNotFoundError, ApplicationError
from skbio.util import create_dir
from bfillings.denoiser import lazy_parse_sff_handle
from burrito.util import which

from qiime.util import (get_qiime_project_dir, FileFormatError)


def write_sff_header(header, fh, num=None):
    """writes the Common header of a sff.txt file.

    header: the header of an sff file as returned by the sff parser.

    fh: output file handle

    num: number of flowgrams to be written in the header.
          Note that his number should match the final number, if
          the resulting sff.txt should be consistent.
          """

    lines = ["Common Header:"]
    if (num is not None):
        header["# of Flows"] = num

    lines.extend(["  %s:\t%s" % (param, header[param])
                  for param in header])
    fh.write("\n".join(lines) + "\n\n")

#    Wrap into explicit function so we can easily move the data dir around.


def get_denoiser_data_dir():
    """Return the directory of the denoiser error profiles.
    """
    dir = get_qiime_project_dir() + "/qiime/support_files/denoiser/Data/"
    return dir


def get_flowgram_ali_exe():
    """Return the executable name of the flowgram alignment prog
    """
    return "FlowgramAli_4frame"


def check_flowgram_ali_exe():
    """Check if we have a working FlowgramAligner"""
    ali_exe = get_flowgram_ali_exe()

    if which(ali_exe) is None:
        raise ApplicationNotFoundError("The alignment program %s is not "
                                       "accessible via the PATH environment "
                                       "variable." % ali_exe)

    # test if its callable and actually works
    command = "%s -h" % ali_exe
    proc = Popen(command, shell=True, universal_newlines=True,
                 stdout=PIPE, stderr=STDOUT)

    if (proc.wait() != 0):
        raise ApplicationError("Calling %s failed. Check permissions and that it is in fact an executable."
                               % ali_exe)

    result = proc.stdout.read()
    # check that the help string looks correct
    if (not result.startswith("Usage")):
        raise ApplicationError("Calling %s failed. Check permissions and that it is in fact an executable."
                               % ali_exe)
    return True


class FlowgramContainerFile():

    """A Flogram container using a file.

    This class can be used to store intermediate flowgram files on disk.
    Slower, but keeps a very low memory footprint.
    """

    def __init__(self, header, outdir="/tmp/"):
        # set up output file
        fd, self.filename = mkstemp(dir=outdir, prefix="fc",
                                   suffix=".sff.txt")
        close(fd)
        self.fh = open(self.filename, "w")
        write_sff_header(header, self.fh)

        self.write_mode = True

    def add(self, flowgram):
        if self.write_mode:
            self.fh.write(flowgram.createFlowHeader() + "\n")
        else:
            raise ValueError(
                "add function can't be called after iteration started.")

    def __iter__(self):
        # make it read_only and reset to start of file
        self.write_mode = False
        self.fh.close()
        (self.flowgrams, self.header) = lazy_parse_sff_handle(
            open(self.filename))
        return self.flowgrams

    def __del__(self):
        remove(self.filename)


class FlowgramContainerArray():

    """A Flogram container using a simple list.

    Keeps all flowgrams in memory. Faster, but needs a lot of memory.
    """

    def __init__(self, header=None):
        """
        header parameter only for compatibility with  FCArray
        """

        self.data = []

    def add(self, flowgram):
        self.data.append(flowgram)

    def __iter__(self):
        return self.data.__iter__()


def make_stats(mapping):
    """Calculates some cluster statistics (counts).

    mapping: The prefix mapping dict
    """
    stats = ["Clustersize\t#"]
    counts = defaultdict(int)
    for key in mapping.keys():
        counts[len(mapping[key])] += 1

    keys = sorted(counts.keys())
    for key in keys:
        stats.append("%d:\t\t%d" % (key + 1, counts[key]))
    return "\n".join(stats)


def sort_ids(ids, mapping):
    """sorts ids based on their cluster_size"""

    def _lookup(id):
        try:
            return len(mapping[id])
        except (KeyError):
            return 0

    deco = [(_lookup(id), id) for id in ids]
    deco.sort(reverse=True)
    return [id for _, id in deco]


def sort_seqs_by_clustersize(seqs, mapping):
    """sort sequences by the cluster size from mapping

    seqs: seqs as iterator or list of (label, seq)

    mapping: cluster mapping as dict
    """
    ids = []
    seqs_cache = {}
    for header, seq in seqs:
        id = header.split("|")[0]
        id = id.rstrip(" ")
        ids.append(id)
        seqs_cache[id] = (header, seq)

    for id in sort_ids(ids, mapping):
        yield seqs_cache[id]


def get_representatives(mapping, seqs):
    """Returns representative seqs.

    mapping: The prefix mapping dict

    seqs_fh: An open Fasta filehandle
    """
    for (label, seq) in seqs:
        if(label in mapping):
            seq = BiologicalSequence(
                seq, id="%s: %d" % (label, len(mapping[label]) + 1))
            yield seq.upper()


def store_mapping(mapping, outdir, prefix):
    """Store the mapping of denoised seq ids to input ids."""
    fh = open(outdir + "/" + prefix + "_mapping.txt", "w")
    for (key, valuelist) in mapping.iteritems():
        fh.write("%s:" % key)
        for v in valuelist:
            fh.write("\t%s" % v)
        fh.write("\n")
    fh.close()


def store_clusters(mapping, sff_fp, outdir="/tmp/", store_members=False):
    """Stores fasta and flogram file for each cluster."""

    # get mapping read to cluster
    invert_map = invert_mapping(mapping)
    (flowgrams, header) = lazy_parse_sff_handle(open(sff_fp))

    leftover_fasta_fh = open(outdir + "/singletons.fasta", "w")
    centroids = []
    for f in flowgrams:
        try:
            key = invert_map[f.Name]
        except KeyError:
            # this flowgram has not been clustered
            continue
        if (len(mapping[key]) == 0):
            # do not store singletons in a separate cluster
            leftover_fasta_fh.write(f.toFasta() + "\n")
            continue
        elif(f.Name in mapping):
            # save as a centroid
            centroids.append((len(mapping[f.Name]) + 1, f.Name, f.toSeq()))

        if (store_members):
            flows_fh = open(outdir + key + ".flows", "a")
            fasta_fh = open(outdir + key + ".fasta", "a")
            flows_fh.write("%s\n" % f)
            fasta_fh.write(f.toFasta() + "\n")
            fasta_fh.close()
            flows_fh.close()

    leftover_fasta_fh.close()

    # sort and store ordered by cluster_size
    centroids.sort(reverse=True)
    centroid_fh = open(outdir + "/centroids.fasta", "w")
    for size, name, seq in centroids:
        centroid_fh.write(">%s | cluster size: %d \n%s\n" %
                          (name, size, seq))
    centroid_fh.close()


def squeeze_seq(seq):
    """Squeezes consecutive identical nucleotides to one.

    seq: a string
    """

    return sub(r'([AGCTacgt])\1+', '\\1', seq)


def wait_for_file(filename, interval=10, test_mode=False):
    """Puts the process to sleep until the file is there.

    filename: file to wait for

    interval: sleep interval in seconds

    test_mode: raise Exception instead of going to sleep
    """
    while(not exists(filename)):
        if test_mode:
            raise RuntimeWarning
        sleep(interval)

def init_flowgram_file(filename=None, n=0, l=400, prefix="/tmp/"):
    """Opens a file in plain flowgram format and writes header information.

    filename: name of output file

    n: number of flowgrams in the file

    l: length of each flowgram in the file

    prefix: directory prefix

    Returns an open filehandle and the file name.
    """

    if (filename is None):
        fd, filename = mkstemp (dir=prefix, suffix=".dat")
        close(fd)

    fh = open(filename, "w")
    fh.write("%d %d\n" % (n, l))
    return (fh, filename)


def append_to_flowgram_file(identifier, flowgram, fh, trim=False):
    """Adds one flowgram to an open plain flowgram file.

    id: identifier of this flowgram

    flowgram: the flowgram itself

    fh: filehandle to write in

    trim: Boolean flag for quality trimming flowgrams
    """

    if trim:
        flowgram = flowgram.getQualityTrimmedFlowgram()

    # store space separated string representation of flowgram
    if (not hasattr(flowgram, "spaced_flowgram")):
        spaced_flowgram_seq = " ".join(map(str, flowgram.flowgram))
        flowgram.spaced_flowgram = spaced_flowgram_seq
    else:
        spaced_flowgram_seq = flowgram.spaced_flowgram

    fh.write("%s %d %s\n" % (identifier, len(flowgram), spaced_flowgram_seq))


def read_signal_probs(file):
    """Read and check the signal probabilty file"""
    f = open(file)
    lines = f.readlines()
    f.close()

    flow_probs = defaultdict(list)
    flow_logs = defaultdict(list)

    for line in lines:
        if line.startswith('#'):
            continue
        for i, num in enumerate(line.strip().split()[2::2]):
            flow_probs[i].append(float(num))
        for i, num in enumerate(line.strip().split()[1::2]):
            flow_logs[i].append(float(num))

    for p in flow_probs:
        s = sum(flow_probs[p])
        flow_probs[p] = [i / s for i in flow_probs[p]]

    return (flow_probs, flow_logs)


def invert_mapping(mapping):
    """Inverts a dictionary mapping.

    Keys are inserted as a special case:
    Ex: {1:(2,3,4)} ==>{1:1, 2:1, 3:1, 4:1}

    Note: This will overwrite an entry if it is redundant.
    """

    invert_map = {}
    for key in mapping.keys():
        invert_map[key] = key
        for id in mapping[key]:
            invert_map[id] = key
    return invert_map


def cat_sff_files(list_of_file_handles):
    """virtually concat several sff files

    list_of_file_handles: list of open filehandles to .sff.txt files

    returns: flowgram generator, header
    """
    # mimicks lazy_parse_sff_handle on multiple files
    # Move to cogent???
    if (list_of_file_handles == []):
        return [], None
    try:
        flowgrams_and_headers = map(
            lazy_parse_sff_handle,
            list_of_file_handles)
    except ValueError:
        raise FileFormatError('Wrong flogram file format. Make sure you pass the sff.txt format ' +
                              'produced by sffinfo. The binary .sff will not work here.')

    flowgram_iterators = [a for a, b in flowgrams_and_headers]
    return chain(*flowgram_iterators), flowgrams_and_headers[0][1]


def files_exist(comma_sep_fps):
    """check if all files in commasep list exist

    comma_sep_fps: list of filenames as comma-separated string
    """

    filenames = comma_sep_fps.split(",")
    for file in filenames:
        if not exists(file):
            return False
    return True


def read_denoiser_mapping(mapping_fh):
    """read the cluster mapping file handle

    mapping_fh: an open file handle to a cluster file.

    Expected format:

    id1: id2 id3
    id4:
    ...
    """
    denoiser_mapping = {}
    for i, line in enumerate(mapping_fh):
        if line == "":
            continue
        centroid, members = line.split(':')
        denoiser_mapping[centroid] = members.split()
    return denoiser_mapping


def write_checkpoint(current_key, ctr, cluster_mapping,
                     ids, bestscores, order, out_fp):
    """write intermediate results to checkpoint file

    current_key: the identifier of the current denoiser round
    ctr: a uniq counter to label the checkpoint
    cluster_mapping: an intermediate cluster mapping as dict
    ids: the dict of active ids
    order:  a list of ids, which defines the order of which flowgrams are clustered
    bestscores: a dict of
    """

    checkpoint_dir = out_fp + "/checkpoints/"
    if (not exists(checkpoint_dir)):
        create_dir(checkpoint_dir)
    out_fp = checkpoint_dir + "/checkpoint%d.pickle" % ctr
    out_fh = open(out_fp, "w")
    pickle.dump(
        (current_key,
         ctr,
         cluster_mapping,
         ids,
         bestscores,
         order),
        out_fh)

    return out_fp


def read_checkpoint(out_fp):
    """Read in information stored in a checkpoint

    out_fp: The path to the checkpoint dir
    """
    pickle_fh = open(out_fp, "r")
    return pickle.load(pickle_fh)


def sort_mapping_by_size(cluster_mapping):
    """Sort the keys of a dict reative to their values length

    cluster_mapping: dict with mapping as list of ids
    """

    return sorted(cluster_mapping.keys(),
                  cmp=lambda a, b: cmp(len(a), len(b)),
                  key=lambda k: cluster_mapping[k], reverse=True)
