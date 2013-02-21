#!/usr/bin/env python
# File created on 30 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

import os
from shutil import rmtree
from glob import glob
from tarfile import open as open_tarfile
from os.path import join, exists, getsize, split, splitext, dirname
from os import makedirs, system, chdir, getcwd
from numpy import array, absolute
from cogent import LoadTree, LoadSeqs
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from cogent.app.util import ApplicationNotFoundError
from qiime.util import get_tmp_filename
from cogent.parse.binary_sff import parse_binary_sff
from qiime.util import (load_qiime_config, 
                        count_seqs,
                        get_qiime_temp_dir,
                        create_dir)
from qiime.parse import (parse_qiime_parameters,
                         parse_distmat_to_dict,
                         parse_distmat,
                         parse_taxa_summary_table)
from biom.parse import parse_biom_table
from qiime.test import (initiate_timeout,
                        disable_timeout,
                        get_test_data_fps)
from qiime.workflow.util import (call_commands_serially,
                                 no_status_updates,
                                 WorkflowError)
from qiime.workflow.downstream import run_beta_diversity_through_plots


# function to stop/start the timeout with longer
# timer for certain tests
def restart_timeout(seconds):
    disable_timeout()
    initiate_timeout(seconds)

class WorkflowTests(TestCase):
    
    def setUp(self):
        """ """
        self.test_data = get_test_data_fps()
        self.files_to_remove = []
        self.dirs_to_remove = []
        
        # Create example output directory
        tmp_dir = get_qiime_temp_dir()
        self.test_out = get_tmp_filename(tmp_dir=tmp_dir,
                                         prefix='core_qiime_analyses_test_',
                                         suffix='',
                                         result_constructor=str)
        self.dirs_to_remove.append(self.test_out)
        create_dir(self.test_out)
        
        self.qiime_config = load_qiime_config()
        self.params = parse_qiime_parameters({})
        
        initiate_timeout(60)
    
    def tearDown(self):
        """ """
        disable_timeout()
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)
        
    def test_unsupported_options_handled_nicely(self):
        """WorkflowError raised on unsupported option """
        self.params['beta_diversity'] = {'blah':"something-broken"}
        self.assertRaises(WorkflowError,
                          run_beta_diversity_through_plots,
                          self.test_data['biom'][0], 
                          self.test_data['map'][0],
                          self.test_out, 
                          call_commands_serially,
                          self.params,
                          self.qiime_config,
                          tree_fp=self.test_data['tree'][0],
                          parallel=False, 
                          status_update_callback=no_status_updates)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

if __name__ == "__main__":
    main()
