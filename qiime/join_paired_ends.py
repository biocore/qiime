#!/usr/bin/env python
# file: join_paired_ends.py

__author__ = "Mike Robeson"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Mike Robeson"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Mike Robeson"
__email__ = "robesonms@ornl.gov"
__status__ = "Development"

from cogent.parse.fastq import MinimalFastqParser
from cogent.app.fastq_join import run_fastqjoin
from cogent.app.seqprep import run_seqprep
from cogent.app.flash import run_flash
from cogent.app.pandaseq import run_pandaseq

join_method_constructors = {}
join_method_names = {'fastq-join':run_fastqjoin,
                     'SeqPrep':run_seqprep,
                     'flash':run_flash,
                     'pandaseq':run_pandaseq}




