#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@gmail.com"

from pyqi.core.command import (
    Command, CommandIn, CommandOut, ParameterCollection)

from qiime.util import subsample_fasta

class SubsampleFasta(Command):
    BriefDescription = "Randomly subsample sequences from a given fasta file"
    LongDescription = ""
    CommandIns = ParameterCollection([
        CommandIn(Name='percent_subsample', DataType=float,
                  Description='Specify the percentage of sequences to '
                  'subsample', Required=True),
        CommandIn(Name='input_fasta', DataType=file,
                  Description='The input FASTA sequences', Required=True)
    ])

    CommandOuts = ParameterCollection([
        CommandOut(Name="subsampled_sequences", DataType=str,
                   Description="Sequences randomly subsampled from the input "
                   "FASTA file")
    ])

    def run(self, **kwargs):
        percent_subsample = kwargs['percent_subsample']
        input_fasta = kwargs['input_fasta']

        if percent_subsample > 1 or percent_subsample <= 0:
            raise ValueError,('percent_subsample must be in range of 0-1')

        return {'subsampled_sequences': subsample_fasta(input_fasta,
                                                        percent_subsample)}

CommandConstructor = SubsampleFasta
