#!/usr/bin/env python
# File created on 14 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import join, abspath, exists
from os import makedirs

from qiime.parallel.wrapper import ParallelWrapper
from qiime.workflow.util import WorkflowLogger


class ParallelMultipleRarefactions(ParallelWrapper):

    def _construct_job_graph(self, input_fp, output_dir, params):
        """Creates the job workflow graph to execute multiple rarefactions on
        parallel

        Parameters
        ----------
        input_fp : str
            Path to the input OTU table
        output_dir : str
            Path to the output directory
        params : dict
            Parameters to use when calling single_rarefaction.py, in the
            form of {param_name: value}
        """
        # Do the parameter parsing
        min_seqs = params['min']
        max_seqs = params['max']
        step = params['step']
        num_reps = params['num_reps']
        input_fp = abspath(input_fp)
        output_dir = abspath(output_dir)
        lineages_included_str = ('--suppress_lineages_included'
                                 if params['suppress_lineages_included']
                                 else '')
        subsample_multinomial_str = ('--subsample_multinomial'
                                     if params['subsample_multinomial']
                                     else '')

        # Create the log file
        self._logger = WorkflowLogger()

        # Create the output directory if it does not exists
        if not exists(output_dir):
            makedirs(output_dir)

        for depth in range(min_seqs, max_seqs + 1, step):
            for rep_num in range(num_reps):
                output_fp = join(output_dir, 'rarefaction_%d_%d.biom'
                                             % (depth, rep_num))
                cmd = ("single_rarefaction.py -i %s -o %s %s %s -d %s"
                       % (input_fp, output_fp, lineages_included_str,
                          subsample_multinomial_str, depth))
                self._job_graph.add_node("%d_%d" % (depth, rep_num),
                                         job=(cmd,), requires_deps=False)
