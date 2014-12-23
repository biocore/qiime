#!/usr/bin/env python

__author__ = "Mike Robeson"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Mike Robeson"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Mike Robeson"
__email__ = "robesonms@ornl.gov"

from bfillings.fastq_join import FastqJoin, join_paired_end_reads_fastqjoin
from bfillings.seqprep import SeqPrep, join_paired_end_reads_seqprep
import os
import gzip

join_method_constructors = {}
join_method_names = {'fastq-join': join_paired_end_reads_fastqjoin,
                     'SeqPrep': join_paired_end_reads_seqprep}



