#!/usr/bin/env python
from os.path import exists, isdir, splitext, join
from os import makedirs, listdir

from skbio.core.distance import SymmetricDistanceMatrix
from skbio.maths.stats.ordination import PCoA

from qiime.format import format_coords

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Antonio Gonzalez Pena",
               "Catherine Lozupone", "Emily TerAvest",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


def pcoa(lines):
    dist_mtx = SymmetricDistanceMatrix.from_file(lines)
    pcoa_obj = PCoA(dist_mtx)
    pcoa_results = pcoa_obj.scores()

    return format_coords(dist_mtx.ids, pcoa_results.species,
                         pcoa_results.eigvals, pcoa_results.perc_expl)


def multiple_file_pcoa(input_dir, output_dir):
    """perform PCoAs on all distance matrices in the input_dir
    """
    if not exists(output_dir):
        makedirs(output_dir)
    file_names = listdir(input_dir)
    file_names = [fname for fname in file_names
                  if not (fname.startswith('.') or isdir(fname))]

    for fname in file_names:
        infile = join(input_dir, fname)
        with open(infile, 'U') as lines:
            pcoa_res_string = pcoa(lines)

        base_fname, ext = splitext(fname)
        outfile = join(output_dir, 'pcoa_%s.txt' % base_fname)
        with open(outfile, 'w') as outf:
            outf.write(pcoa_res_string)
