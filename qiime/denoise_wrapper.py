#!/usr/bin/env python
from __future__ import division

"""Functions supporting the QIIME denoiser"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from os import system, listdir, remove, rmdir
from os.path import exists, split
from re import search
from itertools import chain

from skbio.parse.sequences import parse_fasta
from bfillings.denoiser import lazy_parse_sff_handle
from burrito.util import ApplicationNotFoundError, ApplicationError

from qiime.util import load_qiime_config
from qiime.denoiser.flowgram_clustering import denoise_seqs


def fast_denoiser(
        sff_fps, fasta_fp, tmp_outdir, num_cpus, primer, verbose=True,
        titanium=False):
    """wrapper function calling methods from the Denoiser package."""
    if num_cpus > 1:
        denoise_seqs(sff_fps, fasta_fp, tmp_outdir,
                     primer=primer, cluster=True, num_cpus=num_cpus,
                     verbose=verbose, titanium=titanium)
    else:
        denoise_seqs(sff_fps, fasta_fp, tmp_outdir, primer=primer,
                     verbose=verbose, titanium=titanium)

    # read centroids and singletons
    centroids = parse_fasta(open(tmp_outdir + "/centroids.fasta"))
    singletons = parse_fasta(open(tmp_outdir + "/singletons.fasta"))

    seqs = chain(centroids, singletons)

    # read mapping
    mapping = {}
    cluster_mapping = open(tmp_outdir + "/denoiser_mapping.txt")
    for i, cluster in enumerate(cluster_mapping):
        cluster, members = cluster.split(':')
        members = members.split()
        clust = [cluster]
        clust.extend(members)
        mapping[i] = clust

    return seqs, mapping


def extract_cluster_size(line):
    """ extract the size of a cluster from the header line.

    line is expected to be of this format:
>GCC6FHY01EQVIC | cluster size: 5
    """
    cluster_size = line.split(":")[-1]

    try:
        cluster_size = int(cluster_size)
    except ValueError:
        return 0
    return cluster_size
