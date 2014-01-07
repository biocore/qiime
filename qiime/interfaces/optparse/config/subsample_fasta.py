#!/usr/bin/env python
from __future__ import division

__author__ = "Adam Robbins-Pianka"
__credits__ = ["Adam Robbins-Pianka"]
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"

from pyqi.core.interfaces.optparse import (OptparseUsageExample,
                                           OptparseOption, OptparseResult)
from pyqi.core.command import (make_command_in_collection_lookup_f,
                               make_command_out_collection_lookup_f)
from qiime.commands.subsample_fasta import CommandConstructor

from pyqi.core.interfaces.optparse.input_handler import file_reading_handler
from pyqi.core.interfaces.optparse.output_handler import write_string

# Convenience function for looking up parameters by name.
cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)

usage_examples = [
    OptparseUsageExample(ShortDesc="Example:",
                         LongDesc="Subsample the seqs.fna file, randomly "
                         "select 5% of the sequences:",
                         Ex="%prog -i $PWD/seqs.fna -p 0.05 -o "
                         "$PWD/subsampled_seqs.fna")
]

inputs = [
    OptparseOption(Parameter=cmd_in_lookup('input_fasta'),
                   Type='existing_filepath',
                   Action='store',
                   Handler=file_reading_handler,
                   ShortName='i',
                   Name='input-fasta-fp',
                   Required=True,
                   Help='The input FASTA file path'
                   ),
    OptparseOption(Parameter=None,
                   Type='new_filepath',
                   Action='store',
                   ShortName='o',
                   Name='output-fasta-fp',
                   Required=True,
                   Help='The output FASTA file path'
                   ),
    OptparseOption(Parameter=cmd_in_lookup('percent_subsample'),
                   Type='float',
                   Action='store',
                   Handler=None,
                   ShortName='p',
                   Name='percent-subsample',
                   Required=True,
                   Help='Specify the percentage of sequences to subsample'
                   )
]

outputs = [
    OptparseResult(Parameter=cmd_out_lookup('subsampled_sequences'),
                    Handler=write_string,
                    InputName='output-fasta-fp')
]
