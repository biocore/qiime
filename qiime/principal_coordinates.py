#!/usr/bin/env python
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import PCoA

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Antonio Gonzalez Pena",
               "Catherine Lozupone", "Emily TerAvest",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


def pcoa(lines):
    """Run PCoA on the distance matrix present on lines"""
    # Parse the distance matrix
    dist_mtx = DistanceMatrix.read(lines)
    # Create the PCoA object
    pcoa_obj = PCoA(dist_mtx)
    # Get the PCoA results and return them
    return pcoa_obj.scores()
